"""Simple Materials Interface"""

"""BLOCK: Packages Import"""
# This block handles the import of necessary packages for the script.
# Micropip is used to install packages in a Pyodide environment.
# ASE (Atomic Simulation Environment) is installed as it is required for manipulating atoms, molecules and crystals.
import micropip

await micropip.install("ase")
print("installed ase")


"""BLOCK: Classes and Definitions"""
# Define indices to identify the substrate and layer from the list of input materials.
# These materials are expected to be provided as primitive cells of a 3D material (substrate)
# and a 2D material (layer), which will form the interface.
SUBSTRATE_INDEX = 0
LAYER_INDEX = 1

# Set Miller indices (h, k, l) for the resulting surface
SURFACE_MILLER_H = 1
SURFACE_MILLER_H = 1
SURFACE_MILLER_K = 1
SURFACE_MILLER_L = 1
SURFACE_VACUUM = 5  # defines the vacuum spacing (in Ångströms) added to the SURFACE unit cell.
SURFACE_NUMBER_OF_LAYERS = 3

# Set parameters of the interface
# Matrices of deformation for the surface and the layer
INTERFACE_SURFACE_V_MATRIX = [[1, 0], [0, 1]]
INTERFACE_LAYER_V_MATRIX = [[1, 0], [0, 1]]
INTERFACE_DISTANCE = 2.0


# DO NOT EDIT code below unless you are developing a new transformation.
# The required imports
from ase.build import surface, supercells
from ase.io import read, write
import io
import numpy as np


# The following functions are used to convert between the POSCAR format and the ASE Atoms object.
def ase_poscar_to_atoms(poscar):
    input = io.StringIO(poscar)
    atoms = read(input, format="vasp")
    return atoms


def ase_atoms_to_poscar(atoms):
    output = io.StringIO()
    write(output, atoms, format="vasp")
    content = output.getvalue()
    output.close()
    return content


def expand_matrix_2x2_to_3x3(matrix_2x2):
    matrix_3x3 = np.identity(3)
    matrix_3x3[0:2, 0:2] = matrix_2x2
    return matrix_3x3


# The `MaterialInterface` class encapsulates the creation and manipulation of a material interface.
# It includes methods to create the structure of the interface, calculate strain and distance between layers.
class MaterialInterface:
    def __init__(self, substrate, material, settings=None):
        self.substrate = substrate
        self.material = material
        self.settings = settings
        if settings:
            for key in self.settings.keys():
                if key in settings:
                    self.settings[key].update(settings[key])
        self.structure = self.create_structure()

    def create_structure(self):
        surface = self.settings["surface"]
        interface = self.settings["interface"]

        self.substrate = surface(
            self.substrate,
            (surface["miller:h"], surface["miller:k"], surface["miller:l"]),
            vacuum=surface["vacuum"],
            layers=surface["number_of_layers"],
        )

        surface_v_matrix = expand_matrix_2x2_to_3x3(interface["surface_v:matrix"])
        layer_v_matrix = expand_matrix_2x2_to_3x3(interface["layer_v:matrix"])

        self.substrate = supercells.make_supercell(self.substrate, surface_v_matrix)
        self.substrate.wrap()
        self.material = supercells.make_supercell(self.material, layer_v_matrix)
        self.original_material = self.material.copy()
        self.material.set_cell(self.substrate.get_cell(), scale_atoms=True)
        self.material.wrap()

        z_offset = self.calculate_distance()
        self.material.positions[:, 2] += z_offset
        interface = self.substrate + self.material
        interface.wrap()
        return interface

    def calculate_strain(self, substrate=None, material=None):
        """Calculates strain for the material layer on the substrate"""

        if substrate is None:
            substrate = self.substrate
        if material is None:
            material = self.original_material

        substrate_cell = substrate.get_cell()
        material_cell = material.get_cell()

        a0 = np.linalg.norm(substrate_cell[0])
        b0 = np.linalg.norm(substrate_cell[1])

        a1 = np.linalg.norm(material_cell[0])
        b1 = np.linalg.norm(material_cell[1])

        strain_a = (a1 - a0) / a0
        strain_b = (b1 - b0) / b0

        return (strain_a, strain_b)

    def calculate_distance(self):
        """Calculates distance between the substrate and the material"""
        interface = self.settings["interface"]
        z_max_substrate = max(self.substrate.positions[:, 2])
        z_min_material = min(self.material.positions[:, 2])
        z_offset = z_max_substrate - z_min_material + interface["distance"]

        return z_offset


# The function 'transform' will serve as the main execution point for the transformation.
# It creates a `MaterialInterface` object with input materials and settings,
# calculates the interface structure, and handles output.
def transform():
    settings = {
        "surface": {
            "miller:h": SURFACE_MILLER_H,
            "miller:k": SURFACE_MILLER_K,
            "miller:l": SURFACE_MILLER_L,
            "vacuum": SURFACE_VACUUM,
            "number_of_layers": SURFACE_NUMBER_OF_LAYERS,
        },
        "interface": {
            "surface_v:matrix": INTERFACE_SURFACE_V_MATRIX,
            "layer_v:matrix": INTERFACE_LAYER_V_MATRIX,
            "distance": INTERFACE_DISTANCE,
        },
    }

    materials = globals()["materials_in"]
    substrate_data = materials[SUBSTRATE_INDEX]
    layer_data = materials[LAYER_INDEX]

    substrate = ase_poscar_to_atoms(substrate_data["poscar"])
    layer = ase_poscar_to_atoms(layer_data["poscar"])

    interface = MaterialInterface(substrate, layer, settings)

    print("Interface structure: ", interface.structure)
    print("Strain alongside lattice a:", interface.calculate_strain().a)
    print("Strain alongside lattice b:", interface.calculate_strain().b)

    globals()["materials_out"] = [
        {
            "poscar": ase_atoms_to_poscar(interface.structure),
        }
    ]

    # Return the globals() dictionary to the platform JS environment.
    return globals()


# Execute the main function to preform the material transformation process.
transform()
