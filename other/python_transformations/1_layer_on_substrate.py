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
NOTE: edit the variables above for your specific use case.
"""

# Indices identify the substrate and layer from the list of input materials under `materials_in` in globals().
SUBSTRATE_INDEX = 0
LAYER_INDEX = 1

SETTINGS = {
    "substrate_surface": {
        # Set Miller indices as a tuple for the resulting substrate surface.
        "miller_indices": (1, 1, 1),
        "vacuum": 5,
        "number_of_layers": 3,
        "superlattice_matrix": [
            [1, 0],
            [0, 1]
        ],
    },
    "layer_surface": {
        # Set Miller indices as a tuple for the resulting layer surface.
        "miller_indices": (0, 0, 1),
        "vacuum": 5,
        "number_of_layers": 1,
        "superlattice_matrix": [
            [1, 0],
            [0, 1]
        ],
    },
    "interface": {
        "distance": 3.0,
    },
    # If True the layer cell and basis vectors will be scaled to fit the substrate cell.
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

# The following 3 util functions are used for convenience purposes.
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
    # Create the surface
    surface_atoms = surface(atoms, miller_indices, number_of_layers, vacuum)
    # Expand and apply the superlattice matrix
    expanded_matrix = expand_matrix_2x2_to_3x3(superlattice_matrix)
    supercell_atoms = make_supercell(surface_atoms, expanded_matrix)
    return supercell_atoms


def calculate_strain_matrix(original_layer_cell, scaled_layer_cell):
    # Calculate the original and scaled norms
    cell1 = np.array(original_layer_cell)[0:2,0:2]
    cell2 = np.array(scaled_layer_cell)[0:2,0:2]
    transformation_matrix = np.linalg.inv(cell1) @ cell2
    difference_matrix = transformation_matrix - np.eye(2)
    return difference_matrix


def create_interface(substrate_poscar, layer_poscar, substrate_surface_settings, layer_surface_settings, interface_settings, scale_layer_to_fit=False):
    substrate_atoms = poscar_to_atoms(substrate_poscar)
    layer_atoms = poscar_to_atoms(layer_poscar)
    
    substrate_supercell = create_surface_and_supercell(
        substrate_atoms,
        miller_indices=substrate_surface_settings["miller_indices"],
        number_of_layers=substrate_surface_settings["number_of_layers"],
        vacuum=substrate_surface_settings["vacuum"],
        superlattice_matrix=substrate_surface_settings["superlattice_matrix"]
    )
    
    layer_supercell = create_surface_and_supercell(
        layer_atoms,
        miller_indices=layer_surface_settings["miller_indices"],
        number_of_layers=layer_surface_settings["number_of_layers"],
        vacuum=layer_surface_settings["vacuum"],
        superlattice_matrix=layer_surface_settings["superlattice_matrix"]
    )
    

    # Calculate strain on the layer
    m = calculate_strain_matrix(substrate_supercell.get_cell(), layer_supercell.get_cell())
    m_percent = m * 100 
    percent_str = np.array2string(m_percent, formatter={'float': '{:0.2f}%'.format})
    print("Strain matrix as percentages:\n", percent_str)
    
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
    Creates a `MaterialInterface` based on the input materials, the InterfaceCreator class, and SETTINGS above.
    Makes the resulting material available to the platform JS environment.

    Params:
        (implicit) globals(): The data passed from the platform JS environment.
            - materials_in: The list of input materials. Materials are represented as JS classes. See
                https://github.com/Exabyte-io/made.js/blob/dev/src/material.js#L293 for more information..

    Returns:
        globals(): The globals() dictionary is returned to the platform JS environment.
            - materials_out: The list of output materials.
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
    scale_layer_to_fit=SETTINGS["scale_layer_to_fit"]
    )

    # Make the data available to the platform JS environment.
    globals()["materials_out"] = [
        {
            "poscar": atoms_to_poscar(interface_structure),
        }
    ]
    return globals()


main()
