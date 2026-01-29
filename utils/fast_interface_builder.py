from __future__ import annotations

from typing import List, Optional, Sequence, Union

import numpy as np
from mat3ra.code.array_with_ids import ArrayWithIds
from mat3ra.esse.models.core.abstract.matrix_3x3 import Matrix3x3Schema
from mat3ra.esse.models.materials_category_components.entities.auxiliary.two_dimensional.supercell_matrix_2d import (
    SupercellMatrix2DSchema,
)
from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.interface.simple import InterfaceAnalyzer
from mat3ra.made.tools.analyze.slab import SlabMaterialAnalyzer
from mat3ra.made.tools.build.compound_pristine_structures.two_dimensional.interface.base.build_parameters import (
    InterfaceBuilderParameters,
)
from mat3ra.made.tools.build.compound_pristine_structures.two_dimensional.interface.base.builder import InterfaceBuilder
from mat3ra.made.tools.build.compound_pristine_structures.two_dimensional.interface.base.configuration import (
    InterfaceConfiguration,
)
from mat3ra.made.tools.build_components import MaterialWithBuildMetadata
from mat3ra.made.tools.build_components.entities.core.two_dimensional.vacuum.configuration import VacuumConfiguration


def _embed_2x2_in_3x3(mat_2x2, *, dtype=float) -> np.ndarray:
    m = np.eye(3, dtype=dtype)
    m[:2, :2] = np.array(mat_2x2, dtype=dtype)
    return m


def build_interface_from_2d_match(
    *,
    film_slab: Union[Material, MaterialWithBuildMetadata],
    substrate_slab: Union[Material, MaterialWithBuildMetadata],
    film_supercell_matrix_2x2,
    substrate_supercell_matrix_2x2,
    film_deformation_2x2,
    substrate_deformation_2x2=None,
    align_rotation_2x2=None,
    twist_deg: float = 0.0,
    gap: float = 3.0,
    vacuum: float = 10.0,
    xy_shift: Optional[Sequence[float]] = None,
    reduce_result_cell_to_primitive: bool = False,
) -> Material:
    """
    Build a periodic interface from 2D matching output.

    The match is defined by:
      - Integer in-plane supercell matrices (2x2) for film and substrate
      - In-plane deformation(s) (2x2) applied to film/substrate to make them commensurate
      - Optional alignment rotation (2x2) bringing the substrate first vector to +x
      - Optional in-plane twist (deg) applied to film after alignment

    Returns:
      Mat3ra `Material` interface built via Made's `InterfaceBuilder`.
    """
    xy_shift = list(xy_shift) if xy_shift is not None else [0.0, 0.0]
    if len(xy_shift) != 2:
        raise ValueError("xy_shift must be a 2-element sequence [dx, dy].")

    align_rotation_2x2 = np.array(align_rotation_2x2, dtype=float) if align_rotation_2x2 is not None else np.eye(2)
    substrate_deformation_2x2 = (
        np.array(substrate_deformation_2x2, dtype=float) if substrate_deformation_2x2 is not None else np.eye(2)
    )
    film_deformation_2x2 = np.array(film_deformation_2x2, dtype=float)

    R_align_3 = _embed_2x2_in_3x3(align_rotation_2x2, dtype=float)
    theta = np.deg2rad(float(twist_deg))
    R_twist_3 = _embed_2x2_in_3x3([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]], dtype=float)

    # In Made, xy supercell matrices are 2x2 (applied in-plane only).
    film_xy_supercell = SupercellMatrix2DSchema(root=np.array(film_supercell_matrix_2x2, dtype=int).tolist())
    substrate_xy_supercell = SupercellMatrix2DSchema(root=np.array(substrate_supercell_matrix_2x2, dtype=int).tolist())

    # Strain matrices are 3x3. We fold alignment/twist into the strain so that the final built interface is
    # constructed directly in the matching reference frame.
    #
    # Convention (row-vector lattice matrices): new_cell = old_cell @ strain_matrix
    substrate_strain_3 = R_align_3 @ _embed_2x2_in_3x3(substrate_deformation_2x2, dtype=float)
    film_strain_3 = R_align_3 @ R_twist_3 @ _embed_2x2_in_3x3(film_deformation_2x2, dtype=float)

    substrate_analyzer = SlabMaterialAnalyzer(material=substrate_slab)
    film_analyzer = SlabMaterialAnalyzer(material=film_slab)

    analyzer = InterfaceAnalyzer(
        substrate_slab_configuration=substrate_analyzer.build_configuration,
        film_slab_configuration=film_analyzer.build_configuration,
        substrate_build_parameters=substrate_analyzer.build_parameters,
        film_build_parameters=film_analyzer.build_parameters,
    )

    matched = analyzer.create_matched_configuration_holder(
        substrate_slab_config=substrate_analyzer.build_configuration,
        film_slab_config=film_analyzer.build_configuration,
        match_id=0,
        substrate_xy_supercell_matrix=substrate_xy_supercell,
        film_xy_supercell_matrix=film_xy_supercell,
        substrate_strain_matrix=Matrix3x3Schema(root=substrate_strain_3.tolist()),
        film_strain_matrix=Matrix3x3Schema(root=film_strain_3.tolist()),
        total_strain_percentage=None,
    )

    vacuum_configuration = VacuumConfiguration(size=vacuum)
    config = InterfaceConfiguration(
        stack_components=[matched.substrate_configuration, matched.film_configuration, vacuum_configuration],
        xy_shift=list(xy_shift),
        gaps=ArrayWithIds.from_values([gap]),
    )
    builder = InterfaceBuilder(
        build_parameters=InterfaceBuilderParameters(make_primitive=reduce_result_cell_to_primitive)
    )
    return builder.get_material(config)
