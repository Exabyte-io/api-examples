from typing import Any, Dict, List, Sequence

from mat3ra.made.material import Material


def build_eos_candidate_records(materials: Sequence[Material], scale_factors: Sequence[float]) -> List[Dict[str, Any]]:
    """
    Builds summary records for scaled EOS materials.
    """
    records = []
    for scale_factor, material in zip(scale_factors, materials):
        records.append(
            {
                "scale_factor": scale_factor,
                "material_name": material.name,
                "volume": material.lattice.cell_volume,
            }
        )
    return records
