"""Algorithm from AI"""

import io
from ase.build import surface, make_supercell
from ase.io import read, write
from ase import Atoms
import numpy as np

# Define a function to read POSCAR content and return an Atoms object
def poscar_to_atoms(poscar_content):
    poscar_io = io.StringIO(poscar_content)
    atoms = read(poscar_io, format="vasp")
    return atoms

# Define a function to create a surface and supercell
def prepare_surface_and_supercell(atoms, miller_indices, layers, vacuum, superlattice_matrix):
    # Create the surface
    surface_atoms = surface(atoms, miller_indices, layers, vacuum)
    # Create the supercell
    supercell_atoms = make_supercell(surface_atoms, superlattice_matrix)
    return supercell_atoms

# Define the main function to create the interface
def create_interface(substrate_poscar, layer_poscar, surface_params, superlattice_matrices, z_offset):
    # Convert POSCAR content to Atoms objects
    substrate_atoms = poscar_to_atoms(substrate_poscar)
    layer_atoms = poscar_to_atoms(layer_poscar)
    
    # Prepare the substrate surface and supercell
    substrate_supercell = prepare_surface_and_supercell(
        substrate_atoms, 
        (surface_params["miller"]["h"], surface_params["miller"]["k"], surface_params["miller"]["l"]),
        surface_params["number_of_layers"], 
        surface_params["vacuum"], 
        superlattice_matrices["slab"]
    )
    
    # Prepare the layer supercell
    layer_supercell = prepare_surface_and_supercell(
        layer_atoms, 
        (1, 1, 1),  # Assuming the layer is already a 2D material, just a single layer
        1,  # Only one layer needed for a 2D material
        0,  # No additional vacuum needed, it's already in the POSCAR
        superlattice_matrices["layer"]
    )

    # Adjust Z position based on the Z offset
    z_max_substrate = max(substrate_supercell.positions[:, 2])
    layer_supercell.positions[:, 2] += z_max_substrate + z_offset
    
    # Combine substrate and layer into one Atoms object
    interface = substrate_supercell + layer_supercell
    
    return interface

# Example usage
substrate_poscar = "Your substrate POSCAR content here"
layer_poscar = "Your layer POSCAR content here"
surface_params = {
    "miller": {"h": 1, "k": 1, "l": 1},
    "number_of_layers": 3,
    "vacuum": 10,
}
superlattice_matrices = {
    "slab": [[2, 0, 0], [0, 2, 0], [0, 0, 1]],  # Example superlattice matrix for substrate
    "layer": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],  # Example superlattice matrix for layer
}
z_offset = 3.0  # Example Z offset

# Create the interface
interface_atoms = create_interface(substrate_poscar, layer_poscar, surface_params, superlattice_matrices, z_offset)

# Optional: Save the interface to a file
write("interface_POSCAR.vasp", interface_atoms, format="vasp")

# Print information or visualize the interface
print(interface_atoms)
# view(interface_atoms)  # Uncomment if you have the appropriate visualization tools installed
