import pymatgen.core.surface
import pymatgen.io.ase
import pymatgen.symmetry.analyzer
import ase.build
import ase.constraints

from typing import Tuple, Union, Dict, Any


# Pymatgen

def is_symmetric(slab: pymatgen.core.structure.Structure) -> bool:
    """
    Checks whether a slab is in a point group with inversion symmetry, which includes the following groups:
    -1, 2/m, mmm, 4/m, 4/mmm,-3, -3m, 6/m, 6/mmm, m-3, m-3m
    This technique is borrowed from the `nonstoichiometric_symmetrized_slab` method in PyMatGen's SlabGenerator.

    Args:
        slab (pymatgen.core.structure.Structure): Slab of interest

    Returns:
        True if the slab's spacegroup has inversion symmetry, otherwise False.
    """
    # Create a Spacegroup analyzer object with the slab
    spacegroup = pymatgen.symmetry.analyzer.SpacegroupAnalyzer(slab)
    # Check for inversion symmetry
    is_symmetric = spacegroup.is_laue()
    return is_symmetric


def get_all_slabs_and_terms(crystal: pymatgen.core.structure.Structure, thickness: Union[int, float],
                            is_by_layers: bool) -> Dict[str, Dict[str, Dict[str, Any]]]:
    """
    Gets all slabs and terminations for a given crystal, forcing a specific number of layers in the resultant slab.

    Args:
        crystal (pymatgen.core.structure.Structure):  Crystal of interest
        thickness (int or float): How thick the slab is supposed to be, by either angstroms or number of layers
        is_by_layers (bool): Whether thickness is by the number of layers or by angstroms

    Returns:
        Dict
    """
    # First, get a list of all non-symetrically-equivalent indices that the crystal has
    all_indices = pymatgen.core.surface.get_symmetrically_distinct_miller_indices(crystal, max_index=3)
    slabs = {}

    for plane in all_indices:
        slab_generator = pymatgen.core.surface.SlabGenerator(crystal, plane, min_slab_size=thickness,
                                                             min_vacuum_size=10, center_slab=True,
                                                             in_unit_planes=is_by_layers)

        # Generate all surface terminations, and add them to the list of slabs returned
        all_terminations = slab_generator.get_slabs()
        term_dict = {}
        symmetric_terminations = filter(is_symmetric, all_terminations)
        for term, surface in enumerate(symmetric_terminations):
            ase_surface = pymatgen.io.ase.AseAtomsAdaptor.get_atoms(surface)
            term_dict[str(term)] = {"slab": ase_surface}
        slabs["".join(map(str, plane))] = term_dict
    return slabs


# ASE

def get_bulk_bottom_and_top_frac_coords(slab: ase.Atoms, layers: int = 3) -> Tuple[float, float]:
    """
    Finds the top and bottom of the bulk, in fractional coordinates

    Args:
        slab (ase.Atoms): The slab of interest
        layers (int): How many layers are in the slab

    Returns:
        A tuple containing, in order, the fractional coordintes of the bottom and top of the bulk portion of the slab
    """
    # Determine how large the slab actually is
    c_direction_coords = [atom.scaled_position[2] for atom in slab]
    slab_range = max(c_direction_coords) - min(c_direction_coords)

    # Do some math to figure out where the bottom of the bulk starts
    layer_size = slab_range / layers
    bulk_bottom = min(c_direction_coords) + layer_size

    # Figure out where the top of the bulk starts
    bulk_top = min(c_direction_coords) + 2 * layer_size

    return (bulk_bottom, bulk_top)


def freeze_center_bulk(slab: ase.Atoms) -> None:
    """
    Applies an ASE FixAtoms constraint to atoms found in the center of the slab, in-place.

    Args:
        slab (ase.Atoms): The slab of interest

    Returns:
        None, this function changes the slab in-place.
    """
    # Get the fractional coordinates for the bottom and top of the bulk
    bulk_bottom, bulk_top = get_bulk_bottom_and_top_frac_coords(slab)

    # Filter to get the atoms between the bottom/top of the bulk - that is, the bulk atoms
    frozen_atoms = filter(lambda atom: bulk_bottom <= atom.scaled_position[2] <= bulk_top, slab)

    # Get the indices of the bulk atoms, and apply the constraint to the Atoms object
    frozen_atoms_indices = [atom.index for atom in frozen_atoms]
    fix_atoms_constraint = ase.constraints.FixAtoms(indices=frozen_atoms_indices)
    slab.set_constraint(fix_atoms_constraint)
