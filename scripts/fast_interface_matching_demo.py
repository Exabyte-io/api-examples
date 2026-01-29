from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from utils.interface_matching import FastCoherentMatcher2D


def _hex_vectors(a: float) -> np.ndarray:
    return np.array([[a, 0.0], [0.5 * a, np.sqrt(3.0) * 0.5 * a]])


def _oblique_vectors(a: float, b: float, gamma_deg: float) -> np.ndarray:
    g = np.deg2rad(gamma_deg)
    return np.array([[a, 0.0], [b * np.cos(g), b * np.sin(g)]])


def main() -> None:
    film = _hex_vectors(3.45)
    substrate = _oblique_vectors(4.15, 6.71, 36.0)

    matcher = FastCoherentMatcher2D(
        film,
        substrate,
        max_area=400.0,
        max_results=15,
        # Use relaxed constraints for this "random lattices" example.
        # Tighten these for real interfaces (e.g. 2-5% strain, low shear, small extra rotation).
        max_abs_principal_strain=0.10,
        max_shear_strain=0.10,
        max_extra_rotation_deg=30.0,
        max_length_rel_tol=0.15,
        max_angle_deg_tol=10.0,
        length_bin=0.5,
        angle_bin_deg=2.0,
        shape_weight=0.02,
    )

    for twist in [0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0]:
        matches = matcher.find_matches(twist_deg=twist, strain_share_to_film=1.0, allow_rotation_relaxation=True)
        print(f"\nTwist = {twist:.1f} deg: {len(matches)} matches")
        for m in matches[:5]:
            print(
                "  "
                f"area={m.match_area:8.2f}  "
                f"max|E|={100*m.max_abs_principal_strain:6.2f}%  "
                f"shear={100*m.shear_strain:6.2f}%  "
                f"rot={m.extra_rotation_deg:5.2f}deg  "
                f"Hf={m.film_supercell_matrix.tolist()}  "
                f"Hs={m.substrate_supercell_matrix.tolist()}"
            )


if __name__ == "__main__":
    main()
