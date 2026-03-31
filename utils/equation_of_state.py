from typing import Any, Dict, List, Sequence

from mat3ra.esse.models.core.abstract.matrix_3x3 import Matrix3x3Schema
from mat3ra.made.material import Material
from mat3ra.made.tools.operations.core.unary import strain


def get_material_label(material: Material) -> str:
    """
    Returns a readable label for a material.
    """
    if material.name:
        return material.name
    if getattr(material, "formula", None):
        return material.formula
    return "Material"


def get_isotropic_strain_matrix(scale_factor: float) -> Matrix3x3Schema:
    """
    Returns a 3x3 isotropic strain matrix for uniform lattice scaling.
    """
    return Matrix3x3Schema(
        root=[
            [scale_factor, 0.0, 0.0],
            [0.0, scale_factor, 0.0],
            [0.0, 0.0, scale_factor],
        ]
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

    strain_matrix = get_isotropic_strain_matrix(scale_factor)
    scaled_material = strain(material, strain_matrix)

    base_name = get_material_label(material)
    scaled_material.name = name_template.format(base_name=base_name, scale_factor=scale_factor)
    scaled_material.metadata = dict(material.metadata or {})
    scaled_material.metadata["equationOfState"] = {
        "parentMaterialHash": material.hash,
        "parentMaterialName": base_name,
        "latticeScaleFactor": scale_factor,
        "strainMatrix": strain_matrix.model_dump(),
        "volume": scaled_material.lattice.cell_volume,
    }
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
                "material_name": get_material_label(material),
                "volume": material.lattice.cell_volume,
                "atoms": material.basis.number_of_atoms,
            }
        )
    return records
