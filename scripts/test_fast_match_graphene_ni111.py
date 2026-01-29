from __future__ import annotations

import json
import os
import sys
from math import degrees
from pathlib import Path

import numpy as np

# Ensure local `utils/` is imported (not a site-packages `utils`).
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

# Avoid slow/unwritable matplotlib cache dir warnings triggered by downstream imports.
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

from mat3ra.made.material import Material
from mat3ra.made.tools.build.pristine_structures.two_dimensional.slab import SlabBuilder
from mat3ra.made.tools.build.pristine_structures.two_dimensional.slab.configuration import SlabConfiguration

from utils.interface_matching import FastCoherentMatcher2D, rotation_matrix_2d


def align_first_vector_to_x_right_handed(vectors_2x2: np.ndarray):
    v = np.array(vectors_2x2, dtype=float)
    if np.linalg.det(v) < 0:
        v[1] *= -1
    a = v[0]
    phi = float(np.degrees(np.arctan2(a[1], a[0])))
    R = rotation_matrix_2d(-phi)
    return v @ R, R


def load_material(path: str) -> Material:
    data = json.loads(Path(path).read_text())
    return Material.create(data)


def main() -> None:
    # Input materials shipped with this repo
    substrate = load_material("other/materials_designer/uploads/0-Ni.json")
    film = load_material("other/materials_designer/uploads/1-Graphene.json")

    # Ni(111) substrate slab + graphene(001) film slab
    substrate_slab_config = SlabConfiguration.from_parameters(
        material_or_dict=substrate,
        miller_indices=(1, 1, 1),
        number_of_layers=3,
        vacuum=0.0,
        termination_top_formula=None,
        use_conventional_cell=True,
    )
    film_slab_config = SlabConfiguration.from_parameters(
        material_or_dict=film,
        miller_indices=(0, 0, 1),
        number_of_layers=1,
        vacuum=0.0,
        termination_bottom_formula=None,
        use_conventional_cell=True,
    )

    substrate_slab = SlabBuilder().get_material(substrate_slab_config)
    film_slab = SlabBuilder().get_material(film_slab_config)

    film_vectors = np.array(film_slab.lattice.vector_arrays[0:2])[:, :2]
    substrate_vectors = np.array(substrate_slab.lattice.vector_arrays[0:2])[:, :2]

    substrate_vectors_aligned, R_align = align_first_vector_to_x_right_handed(substrate_vectors)
    film_vectors_aligned = film_vectors @ R_align

    # NOTE:
    # - ZSL's `max_angle_tol` is often treated as radians in practice (e.g. 0.01 rad ~ 0.57 deg).
    # - FastCoherentMatcher2D uses degrees for `max_angle_deg_tol`.
    zsl_like_max_angle_tol = 0.02
    max_angle_deg_tol = max(2.0, degrees(zsl_like_max_angle_tol))

    # Demonstrate the root cause of "0 matches" if you treat graphene's in-plane basis angle as fixed:
    # graphene slab basis here is 120°, while Ni(111) slab basis is 60°. They represent the *same hex lattice*
    # (120° is just an obtuse choice of basis). The matcher must canonicalize that or it will reject on angle.
    strict_matcher = FastCoherentMatcher2D(
        film_vectors_aligned,
        substrate_vectors_aligned,
        canonicalize_obtuse_angle=False,
        max_area=50.0,
        max_results=5,
        max_abs_principal_strain=0.10,
        max_shear_strain=0.10,
        max_extra_rotation_deg=5.0,
        max_length_rel_tol=0.05,
        max_angle_deg_tol=max_angle_deg_tol,
        length_bin=0.25,
        angle_bin_deg=1.0,
    )
    strict_results = strict_matcher.find_matches(
        twist_deg=0.0, strain_share_to_film=1.0, allow_rotation_relaxation=True
    )
    print(f"[strict canonicalize_obtuse_angle=False] matches={len(strict_results)} (expected: 0 for Gr/Ni(111))")

    matcher = FastCoherentMatcher2D(
        film_vectors_aligned,
        substrate_vectors_aligned,
        canonicalize_obtuse_angle=True,
        max_area=50.0,
        max_results=10,
        max_abs_principal_strain=0.05,  # 5%
        max_shear_strain=0.05,  # 5%
        max_extra_rotation_deg=5.0,
        max_length_rel_tol=0.05,
        max_angle_deg_tol=max_angle_deg_tol,
        length_bin=0.25,
        angle_bin_deg=1.0,
    )

    results = matcher.find_matches(twist_deg=0.0, strain_share_to_film=1.0, allow_rotation_relaxation=True)
    print(f"\n[canonicalize_obtuse_angle=True] matches={len(results)} for Gr/Ni(111) at twist=0 deg (max_area=50 Å^2).")
    if not results:
        print("No matches found. Most common causes:")
        print("  - basis angle mismatch (120° vs 60°) if canonicalization is off")
        print("  - overly tight tolerances / too small max_area")
        raise SystemExit(1)

    # Find the expected low-strain commensurate solution:
    # graphene 2x2 on Ni(111) 1x1 (area ~ 21.29 Å^2), with ~0.57% biaxial strain.
    expected = None
    for r in results:
        if (
            r.substrate_supercell_matrix.tolist() == [[1, 0], [0, 1]]
            and abs(np.linalg.det(r.film_supercell_matrix)) == 4
        ):
            expected = r
            break

    best = results[0]
    print("\nTop matches:")
    for r in results[:5]:
        print(
            "  "
            f"area={r.match_area:7.3f}  "
            f"max|E|={100*r.max_abs_principal_strain:6.3f}%  "
            f"Hf={r.film_supercell_matrix.tolist()}  "
            f"Hs={r.substrate_supercell_matrix.tolist()}"
        )

    if expected is None:
        print("\nDid not find graphene 2x2 / Ni(111) 1x1 within constraints.")
        raise SystemExit(2)

    print("\nExpected match (graphene 2x2 on Ni(111) 1x1):")
    print("  film H:", expected.film_supercell_matrix.tolist())
    print("  sub  H:", expected.substrate_supercell_matrix.tolist())
    print(f"  area (Å^2): {expected.match_area:.6f}")
    print(f"  max|principal strain| (%): {100*expected.max_abs_principal_strain:.3f}")
    print(f"  von Mises strain (%): {100*expected.von_mises_strain:.3f}")
    print(f"  extra rotation (deg): {expected.extra_rotation_deg:.6g}")
    print("  deformation (film->substrate):")
    print(expected.film_to_substrate_deformation)


if __name__ == "__main__":
    main()
