from __future__ import annotations

"""
Fast 2D coherent interface (supercell) matching.

This module is intended for interactive workflows where you repeatedly change the in-plane twist angle
and want matches in seconds, with explicit constraints (max strain, max shear, max area, cell shape).

Conventions (match Mat3ra/Made's `np.linalg.solve(film_vectors, substrate_vectors)` usage):
  - In-plane lattice vectors are represented as a 2x2 matrix with row vectors [a; b] in XY.
  - A supercell matrix H (2x2 int) acts on the *left*: `supercell_vectors = H @ vectors`.
  - A Cartesian deformation gradient F (2x2) acts on the *right*: `film_vectors @ F = substrate_vectors`.

Key idea vs. ZSL:
  - Enumerate unique supercells via 2D HNF (det-bounded) deterministically.
  - Use a coarse (length, length, angle) bin index to prune candidate pairs.
  - Compute full deformation/strain only for plausible candidates.

See `scripts/fast_interface_matching_demo.py` for a minimal runnable example.
"""

from dataclasses import dataclass
from functools import cached_property
from math import atan2, ceil, cos, degrees, radians, sin
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np


def _as_2x2_xy(vectors: np.ndarray) -> np.ndarray:
    """
    Convert in-plane lattice vectors to a 2x2 matrix with row vectors (a, b) in XY.

    Accepts:
      - shape (2, 2): already XY
      - shape (2, 3): drops Z
      - shape (3, 3): takes first 2 rows/cols as XY of a,b
    """
    a = np.asarray(vectors, dtype=float)
    if a.shape == (2, 2):
        return a
    if a.shape == (2, 3):
        return a[:, :2]
    if a.shape == (3, 3):
        return a[:2, :2]
    raise ValueError(f"Expected vectors with shape (2,2), (2,3), or (3,3); got {a.shape}")


def rotation_matrix_2d(deg: float) -> np.ndarray:
    theta = radians(deg)
    c = cos(theta)
    s = sin(theta)
    return np.array([[c, -s], [s, c]], dtype=float)


def _cell_area_2d(vectors_2x2: np.ndarray) -> float:
    return float(abs(np.linalg.det(vectors_2x2)))


def _lengths_and_angle_deg(vectors_2x2: np.ndarray) -> Tuple[float, float, float]:
    v1 = vectors_2x2[0]
    v2 = vectors_2x2[1]
    l1 = float(np.linalg.norm(v1))
    l2 = float(np.linalg.norm(v2))
    if l1 == 0 or l2 == 0:
        return l1, l2, 0.0
    cosang = float(np.clip(np.dot(v1, v2) / (l1 * l2), -1.0, 1.0))
    ang = degrees(np.arccos(cosang))
    return l1, l2, ang


def gauss_reduce_basis_2d(vectors_2x2: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Gauss-reduce a 2D lattice basis and return (U, reduced_vectors) with integer unimodular U.

    Conventions:
      - vectors_2x2 is a 2x2 matrix with row vectors [a; b].
      - reduced_vectors = U @ vectors_2x2
      - U is integer with det(U)=±1.
    """
    B = np.array(vectors_2x2, dtype=float)
    U = np.eye(2, dtype=int)

    def swap_rows():
        nonlocal B, U
        B = B[[1, 0], :]
        U = U[[1, 0], :]

    def reduce_b_by_a(mu: int):
        nonlocal B, U
        if mu == 0:
            return
        B[1] = B[1] - float(mu) * B[0]
        U[1] = U[1] - mu * U[0]

    # Ensure non-degenerate
    if abs(np.linalg.det(B)) < 1e-12:
        return U, B

    # Make right-handed
    if np.linalg.det(B) < 0:
        B[1] *= -1
        U[1] *= -1

    # Ensure |a| <= |b|
    if np.linalg.norm(B[0]) > np.linalg.norm(B[1]):
        swap_rows()

    # Gauss reduction loop
    for _ in range(32):
        a = B[0]
        b = B[1]
        aa = float(np.dot(a, a))
        if aa < 1e-18:
            break
        mu = int(np.rint(np.dot(b, a) / aa))
        reduce_b_by_a(mu)

        # Keep |a| <= |b|
        if np.linalg.norm(B[0]) > np.linalg.norm(B[1]):
            swap_rows()

        # Stop if reduced
        a = B[0]
        b = B[1]
        if abs(np.dot(a, b)) <= 0.5 * np.dot(a, a) + 1e-12:
            break

    # Prefer an acute / smaller angle representation (helps matching).
    # Instead of flipping b (which would flip handedness), use an integer basis change b := b + a
    # when dot(a,b) < 0. This keeps det(U) = ±1 and preserves a right-handed basis.
    if np.dot(B[0], B[1]) < 0:
        B[1] = B[1] + B[0]
        U[1] = U[1] + U[0]

    return U, B


def canonicalize_basis_2d_keep_a(vectors_2x2: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Canonicalize a 2D basis while keeping the first vector direction intact.

    This is intentionally lightweight for interactive matching:
      - ensures a right-handed basis (det>0) by flipping b if needed
      - if the basis angle is obtuse (dot(a,b)<0), replaces b := b + a (unimodular)

    Returns (U, canonical_vectors) such that canonical_vectors = U @ vectors_2x2 with integer U.
    """
    B = np.array(vectors_2x2, dtype=float)
    U = np.eye(2, dtype=int)

    if abs(np.linalg.det(B)) < 1e-12:
        return U, B

    # Ensure right-handed
    if np.linalg.det(B) < 0:
        B[1] *= -1
        U[1] *= -1

    # Prefer acute representation without flipping handedness: b := b + a if dot(a,b) < 0
    if np.dot(B[0], B[1]) < 0:
        B[1] = B[1] + B[0]
        U[1] = U[1] + U[0]

    return U, B


def _principal_strains_green_lagrange(F: np.ndarray) -> Tuple[float, float]:
    """
    Principal Green–Lagrange strains from deformation gradient F (2x2).
    E = 1/2 (F^T F - I).
    """
    C = F.T @ F
    E = 0.5 * (C - np.eye(2))
    vals = np.linalg.eigvalsh(E)
    return float(vals[0]), float(vals[1])


def _von_mises_strain_2d_from_green_lagrange(F: np.ndarray) -> float:
    """
    A 2D von Mises-like equivalent strain from Green–Lagrange strain tensor.

    Uses the common plane formula:
      eps_vm = sqrt(e11^2 + e22^2 - e11*e22 + 3*e12^2)
    """
    C = F.T @ F
    E = 0.5 * (C - np.eye(2))
    e11 = float(E[0, 0])
    e22 = float(E[1, 1])
    e12 = float(E[0, 1])
    return float(np.sqrt(e11 * e11 + e22 * e22 - e11 * e22 + 3.0 * e12 * e12))


def polar_decomposition_2x2(F: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Right polar decomposition F = R @ U, where:
      - R is a proper rotation (det(R)=+1 if possible)
      - U is symmetric positive definite stretch
    """
    U_svd, s, Vt = np.linalg.svd(F)
    R = U_svd @ Vt
    # Ensure a proper rotation; adjust the SVD sign convention without corrupting singular values.
    if np.linalg.det(R) < 0:
        Vt = Vt.copy()
        Vt[-1, :] *= -1
        R = U_svd @ Vt
    U = Vt.T @ np.diag(s) @ Vt
    U = 0.5 * (U + U.T)
    return R, U


def spd_matrix_power_2x2(U: np.ndarray, p: float) -> np.ndarray:
    """
    Power of a symmetric positive definite (SPD) 2x2 matrix via eigen-decomposition.
    """
    vals, vecs = np.linalg.eigh(U)
    if np.any(vals <= 0):
        raise ValueError("Expected SPD matrix for power; got non-positive eigenvalue(s).")
    vals_p = np.power(vals, p)
    return vecs @ np.diag(vals_p) @ vecs.T


@dataclass(frozen=True)
class Match2D:
    twist_deg: float
    film_supercell_matrix: np.ndarray  # 2x2 int, left-multiplies film in-plane vectors (rows)
    substrate_supercell_matrix: np.ndarray  # 2x2 int, left-multiplies substrate in-plane vectors (rows)
    match_area: float  # area in Å^2, using the (unstrained) substrate supercell
    film_to_substrate_deformation: np.ndarray  # 2x2, right-multiplies film vectors to get substrate vectors
    film_deformation_applied: np.ndarray  # 2x2, deformation applied to film (depends on strain sharing)
    substrate_deformation_applied: np.ndarray  # 2x2, deformation applied to substrate (depends on strain sharing)
    extra_rotation_deg: float  # rotation angle from polar decomposition of the deformation
    max_abs_principal_strain: float  # Green-Lagrange principal strains
    shear_strain: float  # |E12| from Green-Lagrange
    von_mises_strain: float  # 2D equivalent strain (dimensionless)
    shape_penalty: float  # lower is better


def enumerate_hnf_supercells_2x2(det_max: int) -> List[np.ndarray]:
    """
    Enumerate unique 2D supercells (as integer 2x2 matrices) up to a determinant.

    Uses the (column) HNF enumeration for uniqueness, then transposes to match our
    "row-vector" convention where supercell_vectors = H @ vectors.
    """
    if det_max < 1:
        return []

    mats: List[np.ndarray] = []
    for det in range(1, det_max + 1):
        for a in range(1, det + 1):
            if det % a != 0:
                continue
            c = det // a
            # Column-HNF: [[a, b], [0, c]] with 0 <= b < c
            for b in range(0, c):
                h_col = np.array([[a, b], [0, c]], dtype=int)
                mats.append(h_col.T.copy())
    return mats


class FastCoherentMatcher2D:
    """
    Fast, deterministic 2D coherent supercell matcher for interactive use.

    Key choices vs. ZSL:
      - Deterministic HNF supercell enumeration (no flakiness).
      - Explicit, user-controllable constraints on strain/shear/rotation/area/shape.
      - Designed to be called repeatedly for different twist angles.
    """

    def __init__(
        self,
        film_vectors: np.ndarray,
        substrate_vectors: np.ndarray,
        *,
        canonicalize_obtuse_angle: bool = True,
        reduce_input_bases: bool = False,
        max_area: float = 400.0,
        max_results: int = 25,
        max_length_rel_tol: float = 0.05,
        max_angle_deg_tol: float = 2.0,
        max_extra_rotation_deg: float = 1.0,
        max_abs_principal_strain: float = 0.03,
        max_shear_strain: float = 0.02,
        shape_target_angle_deg: float = 90.0,
        shape_target_aspect: float = 1.0,
        shape_weight: float = 0.05,
        length_bin: float = 0.25,
        angle_bin_deg: float = 1.0,
    ):
        film_vectors_2x2 = _as_2x2_xy(film_vectors)
        substrate_vectors_2x2 = _as_2x2_xy(substrate_vectors)

        self.reduce_input_bases = bool(reduce_input_bases)
        self.canonicalize_obtuse_angle = bool(canonicalize_obtuse_angle)

        U_f = np.eye(2, dtype=int)
        U_s = np.eye(2, dtype=int)
        B_f = film_vectors_2x2
        B_s = substrate_vectors_2x2

        if self.canonicalize_obtuse_angle:
            U_cf, B_f = canonicalize_basis_2d_keep_a(B_f)
            U_cs, B_s = canonicalize_basis_2d_keep_a(B_s)
            U_f = U_cf @ U_f
            U_s = U_cs @ U_s

        if self.reduce_input_bases:
            U_rf, B_f = gauss_reduce_basis_2d(B_f)
            U_rs, B_s = gauss_reduce_basis_2d(B_s)
            U_f = U_rf @ U_f
            U_s = U_rs @ U_s

        self._film_basis_transform = U_f
        self._substrate_basis_transform = U_s
        self.film_vectors = B_f
        self.substrate_vectors = B_s

        self.max_area = float(max_area)
        self.max_results = int(max_results)

        self.max_length_rel_tol = float(max_length_rel_tol)
        self.max_angle_deg_tol = float(max_angle_deg_tol)
        self.max_extra_rotation_deg = float(max_extra_rotation_deg)
        self.max_abs_principal_strain_limit = float(max_abs_principal_strain)
        self.max_shear_strain_limit = float(max_shear_strain)

        self.shape_target_angle_deg = float(shape_target_angle_deg)
        self.shape_target_aspect = float(shape_target_aspect)
        self.shape_weight = float(shape_weight)

        self.length_bin = float(length_bin)
        self.angle_bin_deg = float(angle_bin_deg)

        film_area0 = _cell_area_2d(self.film_vectors)
        sub_area0 = _cell_area_2d(self.substrate_vectors)
        if film_area0 <= 0 or sub_area0 <= 0:
            raise ValueError("Film/substrate in-plane vectors must have non-zero area.")

        self._film_det_max = int(self.max_area // film_area0)
        self._sub_det_max = int(self.max_area // sub_area0)

    @cached_property
    def _film_candidates(
        self,
    ) -> List[Tuple[np.ndarray, np.ndarray, np.ndarray, float, Tuple[float, float, float]]]:
        """
        Cached film supercells independent of twist.

        Each entry is:
          (Hf, Bf0, inv(Bf0), det(Bf0), (l1, l2, angle_deg))

        For an in-plane twist R, the twisted basis is Bf = Bf0 @ R and inv(Bf) = R^T @ inv(Bf0).
        """
        mats = enumerate_hnf_supercells_2x2(self._film_det_max)
        out = []
        for hf in mats:
            B0 = hf @ self.film_vectors
            det0 = float(np.linalg.det(B0))
            if abs(det0) < 1e-12:
                continue
            if abs(det0) > self.max_area + 1e-9:
                continue
            try:
                invB0 = np.linalg.inv(B0)
            except np.linalg.LinAlgError:
                continue
            l1, l2, ang = _lengths_and_angle_deg(B0)
            out.append((hf, B0, invB0, det0, (l1, l2, ang)))

        out.sort(key=lambda x: (_cell_area_2d(x[1]), self._shape_penalty(x[1])))
        return out

    @cached_property
    def _substrate_candidates(self) -> List[Tuple[np.ndarray, np.ndarray, float, Tuple[float, float, float]]]:
        mats = enumerate_hnf_supercells_2x2(self._sub_det_max)
        out = []
        for hs in mats:
            Bs = hs @ self.substrate_vectors
            dets = float(np.linalg.det(Bs))
            if abs(dets) < 1e-12:
                continue
            area = abs(dets)
            if area > self.max_area + 1e-9:
                continue
            l1, l2, ang = _lengths_and_angle_deg(Bs)
            out.append((hs, Bs, dets, (l1, l2, ang)))

        out.sort(key=lambda x: (_cell_area_2d(x[1]), self._shape_penalty(x[1])))
        return out

    @cached_property
    def _substrate_bin_index(
        self,
    ) -> Dict[Tuple[int, int, int], List[int]]:
        index: Dict[Tuple[int, int, int], List[int]] = {}
        for i, (_hs, _Bs, _det, (l1, l2, ang)) in enumerate(self._substrate_candidates):
            key = (
                int(round(l1 / self.length_bin)),
                int(round(l2 / self.length_bin)),
                int(round(ang / self.angle_bin_deg)),
            )
            index.setdefault(key, []).append(i)
        return index

    def _shape_penalty(self, vectors_2x2: np.ndarray) -> float:
        l1, l2, ang = _lengths_and_angle_deg(vectors_2x2)
        if l1 == 0 or l2 == 0:
            return float("inf")
        aspect = l1 / l2
        aspect_pen = abs(aspect / self.shape_target_aspect - 1.0)
        angle_pen = abs(ang - self.shape_target_angle_deg) / max(1e-9, self.shape_target_angle_deg)
        return aspect_pen + angle_pen

    def _candidate_substrate_indices_for(self, film_desc: Tuple[float, float, float]) -> Iterable[int]:
        l1f, l2f, angf = film_desc
        dl1 = max(self.max_length_rel_tol * l1f, self.length_bin)
        dl2 = max(self.max_length_rel_tol * l2f, self.length_bin)
        dang = max(self.max_angle_deg_tol, self.angle_bin_deg)

        l1_bins = int(ceil(dl1 / self.length_bin))
        l2_bins = int(ceil(dl2 / self.length_bin))
        a_bins = int(ceil(dang / self.angle_bin_deg))

        k1 = int(round(l1f / self.length_bin))
        k2 = int(round(l2f / self.length_bin))
        ka = int(round(angf / self.angle_bin_deg))

        seen: set[int] = set()
        for i1 in range(k1 - l1_bins, k1 + l1_bins + 1):
            for i2 in range(k2 - l2_bins, k2 + l2_bins + 1):
                for ia in range(ka - a_bins, ka + a_bins + 1):
                    for idx in self._substrate_bin_index.get((i1, i2, ia), []):
                        if idx not in seen:
                            seen.add(idx)
                            yield idx

    def find_matches(
        self,
        *,
        twist_deg: float = 0.0,
        strain_share_to_film: float = 1.0,
        allow_rotation_relaxation: bool = True,
    ) -> List[Match2D]:
        """
        Args:
            twist_deg: In-plane rotation applied to film basis before matching.
            strain_share_to_film: 1.0 => film takes all stretch; 0.0 => substrate takes all stretch.
            allow_rotation_relaxation: if True, uses polar decomposition and does not count the rotation part as strain.
        """
        twist_deg = float(twist_deg)
        share = float(strain_share_to_film)
        if not (0.0 <= share <= 1.0):
            raise ValueError("strain_share_to_film must be within [0, 1].")

        R_twist = rotation_matrix_2d(twist_deg)
        R_twist_T = R_twist.T
        matches: List[Match2D] = []

        for hf, Bf0, invBf0, detBf0, (l1f, l2f, angf) in self._film_candidates:
            Bf = Bf0 @ R_twist
            invBf = R_twist_T @ invBf0
            for idx in self._candidate_substrate_indices_for((l1f, l2f, angf)):
                hs, Bs, detBs, _desc_s = self._substrate_candidates[idx]

                # Try both (a,b) and (b,a) choices for the substrate basis to avoid missing
                # valid matches due to vector ordering.
                for swap in (False, True):
                    if swap:
                        P = np.array([[0, 1], [1, 0]], dtype=int)
                        hs_eff = P @ hs
                        Bs_eff = P @ Bs
                    else:
                        hs_eff = hs
                        Bs_eff = Bs

                    # Keep orientation consistent to avoid artificial rotation/strain.
                    if detBf0 * detBs < 0:
                        hs_eff = hs_eff.copy()
                        Bs_eff = Bs_eff.copy()
                        hs_eff[1, :] *= -1
                        Bs_eff[1, :] *= -1

                    # Deformation that maps film basis -> substrate basis in Cartesian (row-vector convention):
                    #   Bf @ F = Bs  =>  F = inv(Bf) @ Bs
                    F = invBf @ Bs_eff

                    # Optional: ignore pure rotation as "strain" and constrain extra rotation separately.
                    if allow_rotation_relaxation:
                        R_extra, U = polar_decomposition_2x2(F)
                        rot_angle = abs(degrees(atan2(R_extra[1, 0], R_extra[0, 0])))
                        if rot_angle > self.max_extra_rotation_deg:
                            continue

                        # Split only the stretch between film/substrate
                        try:
                            U_f = spd_matrix_power_2x2(U, share)
                            U_s = spd_matrix_power_2x2(U, share - 1.0)
                        except ValueError:
                            continue

                        F_film = R_extra @ U_f
                        F_sub = U_s
                        # Evaluate strain on both parties; accept if both satisfy.
                        e1f, e2f = _principal_strains_green_lagrange(F_film)
                        e1s, e2s = _principal_strains_green_lagrange(F_sub)
                        max_abs_principal = max(abs(e1f), abs(e2f), abs(e1s), abs(e2s))
                        if max_abs_principal > self.max_abs_principal_strain_limit:
                            continue

                        Ef = 0.5 * ((F_film.T @ F_film) - np.eye(2))
                        Es = 0.5 * ((F_sub.T @ F_sub) - np.eye(2))
                        shear = max(abs(float(Ef[0, 1])), abs(float(Es[0, 1])))
                        if shear > self.max_shear_strain_limit:
                            continue

                        vm = max(
                            _von_mises_strain_2d_from_green_lagrange(F_film),
                            _von_mises_strain_2d_from_green_lagrange(F_sub),
                        )
                        shape_pen = self._shape_penalty(Bs_eff)
                        _score = max_abs_principal + self.shape_weight * shape_pen + 1e-6 * _cell_area_2d(Bs_eff)

                        matches.append(
                            Match2D(
                                twist_deg=twist_deg,
                                film_supercell_matrix=(hf @ self._film_basis_transform).copy(),
                                substrate_supercell_matrix=(hs_eff @ self._substrate_basis_transform).copy(),
                                match_area=_cell_area_2d(Bs_eff),
                                film_to_substrate_deformation=F.copy(),
                                film_deformation_applied=F_film.copy(),
                                substrate_deformation_applied=F_sub.copy(),
                                extra_rotation_deg=rot_angle,
                                max_abs_principal_strain=max_abs_principal,
                                shear_strain=shear,
                                von_mises_strain=vm,
                                shape_penalty=shape_pen,
                            )
                        )
                    else:
                        e1, e2 = _principal_strains_green_lagrange(F)
                        max_abs_principal = max(abs(e1), abs(e2))
                        if max_abs_principal > self.max_abs_principal_strain_limit:
                            continue

                        E = 0.5 * ((F.T @ F) - np.eye(2))
                        shear = abs(float(E[0, 1]))
                        if shear > self.max_shear_strain_limit:
                            continue

                        vm = _von_mises_strain_2d_from_green_lagrange(F)
                        shape_pen = self._shape_penalty(Bs_eff)
                        matches.append(
                            Match2D(
                                twist_deg=twist_deg,
                                film_supercell_matrix=(hf @ self._film_basis_transform).copy(),
                                substrate_supercell_matrix=(hs_eff @ self._substrate_basis_transform).copy(),
                                match_area=_cell_area_2d(Bs_eff),
                                film_to_substrate_deformation=F.copy(),
                                film_deformation_applied=F.copy(),
                                substrate_deformation_applied=np.eye(2, dtype=float),
                                extra_rotation_deg=0.0,
                                max_abs_principal_strain=max_abs_principal,
                                shear_strain=shear,
                                von_mises_strain=vm,
                                shape_penalty=shape_pen,
                            )
                        )

        matches.sort(key=lambda m: (m.max_abs_principal_strain, m.shear_strain, m.shape_penalty, m.match_area))
        return matches[: self.max_results]

    def find_matches_over_twist(
        self,
        twist_deg_values: Iterable[float],
        *,
        strain_share_to_film: float = 1.0,
        allow_rotation_relaxation: bool = True,
        per_twist_max_results: int = 5,
        max_total_results: Optional[int] = None,
    ) -> List[Match2D]:
        """
        Convenience helper for the "let twist be anything" workflow.

        Evaluates a set of twist angles and returns the globally best matches across all angles.
        """
        max_total = int(max_total_results) if max_total_results is not None else self.max_results
        per_twist_max = int(per_twist_max_results)
        if per_twist_max < 1:
            raise ValueError("per_twist_max_results must be >= 1.")
        if max_total < 1:
            return []

        original_max_results = self.max_results
        try:
            # Temporarily cap results per twist to keep the sweep fast and predictable.
            self.max_results = per_twist_max
            all_matches: List[Match2D] = []
            for t in twist_deg_values:
                all_matches.extend(
                    self.find_matches(
                        twist_deg=float(t),
                        strain_share_to_film=strain_share_to_film,
                        allow_rotation_relaxation=allow_rotation_relaxation,
                    )
                )
        finally:
            self.max_results = original_max_results

        all_matches.sort(key=lambda m: (m.max_abs_principal_strain, m.shear_strain, m.shape_penalty, m.match_area))
        return all_matches[:max_total]
