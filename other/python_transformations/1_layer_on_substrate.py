"""TITLE: Place a 2D material Layer on a Surface."""

"""BLOCK: Install Dependencies"""
"""
This scripts creates an interface between two materials: (1) substrate and (2) layer.
It is assumed that the substrate is a 3D material and the layer is a 2D material.
This script is based on the ASE package (Atomic Simulation Environment) and the tutorial at:
https://wiki.fysik.dtu.dk/ase/gettingstarted/manipulating_atoms/manipulating_atoms.html#interface-building.
"""

# This block handles the installation of necessary packages for the script using Micropip in a Pyodide environment.
# ASE (Atomic Simulation Environment) is installed for manipulating structures.
import micropip

await micropip.install("ase")
print("installed ase")

"""BLOCK: Utils, Class Definitions, and main()"""
"""
NOTE: edit the variables below for your specific use case.
"""

# Indices identify the substrate and layer from the list of input materials under `materials_in` in globals().
SUBSTRATE_INDEX = 0
LAYER_INDEX = 1

SETTINGS = {
    "substrate_surface": {
        # Set Miller indices as a tuple for the resulting substrate surface.
        "miller_indices": (1, 1, 1),
        # The vacuum space (in Ångströms) added to the surface in the direction perpendicular to the surface.
        "vacuum": 5,
        # The number of atomic layers in the resulting substrate.
        "number_of_layers": 3,
        # The transformation matrix for the surface. Format is: [[v1x, v1y], [v2x, v2y]].
        # fmt: off
        "superlattice_matrix": [
            [1, 0],
            [0, 1]
        ],
        # fmt: on
    },
    "layer_surface": {
        # Set Miller indices as a tuple for the resulting layer surface: (0,0,1) for 2D material
        "miller_indices": (0, 0, 1),
        # The vacuum space (in Ångströms) added to the surface in the direction perpendicular to the surface.
        "vacuum": 5,
        # The number of atomic layers in the resulting substrate: 1 for 2D material
        "number_of_layers": 1,
        # The transformation matrix for the surface. Format is: [[v1x, v1y], [v2x, v2y]].
        # fmt: off
        "superlattice_matrix": [
            [1, 0],
            [0, 1]
        ],
        # fmt: on
    },
    "interface": {
        "distance": 3.0,
    },
    # If True the layer cell and basis vectors will be scaled to fit the substrate cell.
    # Mind the strain that is introduced by this operation.
    "scale_layer_to_fit": False,
}


"""
NOTE: DO NOT edit code below unless you know what you are doing.
"""

from ase.build import surface as make_surface, supercells
from ase.io import read, write
import io
import numpy as np
import io
from ase.build import surface, make_supercell
from ase.io import read, write
import numpy as np


# Utility functions
def poscar_to_atoms(poscar):
    input = io.StringIO(poscar)
    atoms = read(input, format="vasp")

    return atoms


def atoms_to_poscar(atoms):
    output = io.StringIO()
    write(output, atoms, format="vasp")
    content = output.getvalue()
    output.close()

    return content


def expand_matrix_2x2_to_3x3(matrix_2x2):
    matrix_3x3 = np.identity(3)
    matrix_3x3[0:2, 0:2] = matrix_2x2

    return matrix_3x3


def create_surface_and_supercell(atoms, miller_indices, number_of_layers, vacuum, superlattice_matrix):
    """
    Creates a surface and supercell based on the input parameters.

    Params:
        atoms: The input atoms.
        miller_indices: The Miller indices for the surface.
        number_of_layers: The number of layers in the resulting surface.
        vacuum: The vacuum space (in Ångströms) added to the surface in the direction perpendicular to the surface.
        superlattice_matrix: The transformation matrix for the surface. Format is: [[v1x, v1y], [v2x, v2y]].

    Returns:
        supercell_atoms: The resulting supercell.
    """

    surface_atoms = surface(atoms, miller_indices, number_of_layers, vacuum)
    expanded_matrix = expand_matrix_2x2_to_3x3(superlattice_matrix)
    supercell_atoms = make_supercell(surface_atoms, expanded_matrix)

    return supercell_atoms


def calculate_strain_matrix(scaled_layer_cell, original_layer_cell):
    """
    Calculates the strain matrix based on the scaled and original layer cells.

    Params:
        scaled_layer_cell: The scaled layer cell.
        original_layer_cell: The original layer cell.

    Returns:
        difference_matrix: The difference matrix.
    """

    cell1 = np.array(original_layer_cell)[0:2, 0:2]
    cell2 = np.array(scaled_layer_cell)[0:2, 0:2]
    transformation_matrix = np.linalg.inv(cell1) @ cell2
    difference_matrix = transformation_matrix - np.eye(2)

    return difference_matrix


def create_interface(
    substrate_poscar,
    layer_poscar,
    substrate_surface_settings,
    layer_surface_settings,
    interface_settings,
    scale_layer_to_fit=False,
):
    """
    Creates an interface between the substrate and the layer.

    Params:
        substrate_poscar: The POSCAR data for the substrate.
        layer_poscar: The POSCAR data for the layer.
        substrate_surface_settings: The settings for the substrate surface.
        layer_surface_settings: The settings for the layer surface.
        interface_settings: The settings for the interface.
        scale_layer_to_fit: If True the layer cell and basis vectors will be scaled to fit the substrate cell.

    Returns:
        interface: The resulting interface.
    """

    substrate_atoms = poscar_to_atoms(substrate_poscar)
    layer_atoms = poscar_to_atoms(layer_poscar)

    substrate_supercell = create_surface_and_supercell(
        substrate_atoms,
        miller_indices=substrate_surface_settings["miller_indices"],
        number_of_layers=substrate_surface_settings["number_of_layers"],
        vacuum=substrate_surface_settings["vacuum"],
        superlattice_matrix=substrate_surface_settings["superlattice_matrix"],
    )

    layer_supercell = create_surface_and_supercell(
        layer_atoms,
        miller_indices=layer_surface_settings["miller_indices"],
        number_of_layers=layer_surface_settings["number_of_layers"],
        vacuum=layer_surface_settings["vacuum"],
        superlattice_matrix=layer_surface_settings["superlattice_matrix"],
    )

    # Calculate strain on the layer
    m = calculate_strain_matrix(substrate_supercell.get_cell(), layer_supercell.get_cell())
    m_percent = m * 100
    percent_str = np.array2string(m_percent, formatter={"float": "{:0.2f}%".format})
    print("Strain matrix as percentages:\n", percent_str)

    # Scale the layer to fit the substrate (mind the strain)
    if scale_layer_to_fit:
        layer_supercell.set_cell(substrate_supercell.get_cell(), scale_atoms=True)
        layer_supercell.wrap()

    # Adjust Z position based on the Z offset
    z_max_substrate = max(substrate_supercell.positions[:, 2])
    z_min_layer = min(layer_supercell.positions[:, 2])
    z_offset = z_max_substrate - z_min_layer + interface_settings["distance"]
    layer_supercell.positions[:, 2] += z_offset

    # Combine substrate and layer into one Atoms object
    interface = substrate_supercell + layer_supercell
    return interface


def main():
    """
    The main function of the script.
    Creates an interface between the substrate and the layer based on the SETTINGS defined above.
    Makes the resulting material available to the platform JS environment.

    Params:
        (implicit) globals(): The data passed from the platform JS environment.
            - materials_in: The list of input materials. Materials are represented as JS classes. See
                https://github.com/Exabyte-io/made.js/blob/dev/src/material.js#L293 for more information..

    Returns:
        globals(): The globals() dictionary is returned to the platform JS environment.
            - materials_out: The list of output materials with poscar as strucutre representation.
    """
    # Get the input materials from the platform JS environment.
    materials = globals()["materials_in"]

    # Get poscar data for the substrate and the layer.
    # Per https://github.com/Exabyte-io/made.js/blob/dev/src/material.js#L293
    substrate_data = materials[SUBSTRATE_INDEX].getAsPOSCAR()
    layer_data = materials[LAYER_INDEX].getAsPOSCAR()

    # Interface is created based on SETTINGS for both substrate and layer
    interface_structure = create_interface(
        substrate_poscar=substrate_data,
        layer_poscar=layer_data,
        substrate_surface_settings=SETTINGS["substrate_surface"],
        layer_surface_settings=SETTINGS["layer_surface"],
        interface_settings=SETTINGS["interface"],
        scale_layer_to_fit=SETTINGS["scale_layer_to_fit"],
    )

    # Make the data available to the platform JS environment.
    globals()["materials_out"] = [
        {
            "poscar": atoms_to_poscar(interface_structure),
        }
    ]
    return globals()


main()
