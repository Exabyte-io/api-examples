from typing import Any, Dict, List, Sequence

from mat3ra.made.material import Material
from mat3ra.made.tools.build_components.operations.core.modifications.strain.helpers import (
    create_strain,
    get_isotropic_strain_matrix,
)


def scale_material(
    material: Material,
    scale_factor: float,
    name_template: str = "{base_name} (scale={scale_factor:.4f})",
) -> Material:
    """
    Returns a new material with lattice vectors scaled isotropically.

    Args:
        material (Material): Source material.
        scale_factor (float): Multiplicative scale applied to lattice vectors.
        name_template (str): Template for the new material name.

    Returns:
        Material: Scaled material.
    """
    if scale_factor <= 0:
        raise ValueError("scale_factor must be positive.")

    scaled_material = create_strain(material, get_isotropic_strain_matrix(scale_factor))

    base_name = material.name
    scaled_material.name = name_template.format(base_name=base_name, scale_factor=scale_factor)
    return scaled_material


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
