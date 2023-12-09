"""
TITLE: Place a 2D material Layer on a Surface
This scripts creates an interface between two materials: (1) substrate and (2) layer.
It is assumed that the substrate is a 3D material and the layer is a 2D material.
This script is based on the ASE package (Atomic Simulation Environment) and the tutorial at:
https://wiki.fysik.dtu.dk/ase/gettingstarted/manipulating_atoms/manipulating_atoms.html#interface-building.
"""

"""BLOCK: Install Dependencies"""
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
    "surface": {
        # Set Miller indices (h, k, l) for the resulting substrate surface.
        "miller:h": 1,
        "miller:k": 1,
        "miller:l": 1,
        # The vacuum space (in Ångströms) added to the surface in the direction perpendicular to the surface.
        "vacuum": 5,
        # The number of atomic layers in the resulting substrate.
        "number_of_layers": 3,
    },
    "interface": {
        # The transformation matrix for the surface. Format is: [[v1x, v1y], [v2x, v2y]].
        "surface_v:matrix": [
            [1, 0],
            [0, 1]
        ],
        # The transformation matrix for the layer. Format is the same as above.
        "layer_v:matrix": [
            [1, 0],
            [0, 1]
        ],
        # Distance between the substrate and the layer (in Ångströms).
        "distance": 3.0,
    },
}


"""
NOTE: DO NOT edit code below unless you know what you are doing.
"""

from ase.build import surface as make_surface, supercells
from ase.io import read, write
import io
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


class InterfaceCreator:
    """
    This class Encapsulates the creation and manipulation of a material interface.
    Includes methods to create the structure of the interface, calculate strain and distance between layers.
    """
    def __init__(self, substrate, layer, settings=None):
        self.substrate = substrate
        self.layer = layer
        self.original_layer = self.layer.copy()
        self.settings = settings
        if settings:
            for key in self.settings.keys():
                if key in settings:
                    self.settings[key].update(settings[key])

    def create(self):
        """
        Creates the interface structure from the substrate and the material.
        """
        surface = self.settings["surface"]
        interface = self.settings["interface"]

        # First, create the substrate surface.
        self.substrate = make_surface(
            self.substrate,
            (surface["miller:h"], surface["miller:k"], surface["miller:l"]),
            vacuum=surface["vacuum"],
            layers=surface["number_of_layers"],
        )

        surface_v_matrix = expand_matrix_2x2_to_3x3(interface["surface_v:matrix"])
        layer_v_matrix = expand_matrix_2x2_to_3x3(interface["layer_v:matrix"])

        # Then, create the supercells for both substrate and layer.
        self.substrate = supercells.make_supercell(self.substrate, surface_v_matrix)
        self.substrate.wrap()
        self.layer = supercells.make_supercell(self.layer, layer_v_matrix)

        # Scale the layer cell to match the substrate cell.
        self.layer.set_cell(self.substrate.get_cell(), scale_atoms=True)

        # Workaround: the y-axis of the layer is flipped to match the substrate.
        cell = self.layer.get_cell_lengths_and_angles()
        # cell[5] is the angle(a,b) per:
        # https://wiki.fysik.dtu.dk/ase/ase/atoms.html#ase.Atoms.get_cell_lengths_and_angles
        cell[5] = 180 - cell[5]
        self.layer.set_cell(cell, scale_atoms=True)
        self.layer.wrap()

        # At this point, both self.layer and self.substrate have the same cell.
        # Apply the offset for the layer coordinates to place it on top of the substrate.
        z_offset = self.calculate_layer_offset()
        self.layer.positions[:, 2] += z_offset

        # Stack the layer on top of the substrate.
        interface = self.substrate + self.layer
        interface.wrap()

        return interface

    def calculate_strain(self, substrate=None, material=None):
        """
        Calculates strain between the layer and the substrate.
        """
        if substrate is None:
            substrate = self.substrate
        if material is None:
            material = self.original_layer

        substrate_cell = substrate.get_cell()
        material_cell = material.get_cell()

        a0 = np.linalg.norm(substrate_cell[0])
        b0 = np.linalg.norm(substrate_cell[1])

        a1 = np.linalg.norm(material_cell[0])
        b1 = np.linalg.norm(material_cell[1])

        strain_a = (a1 - a0) / a0
        strain_b = (b1 - b0) / b0

        return {"a": strain_a, "b": strain_b}

    def calculate_layer_offset(self):
        """
        Calculates the offset for the layer coordinates to place it on top of the substrate.
        Uses the distance between the layer and the substrate from the settings.
        Assumes that both the layer and the substrate have the same cell.
        """
        interface = self.settings["interface"]
        z_max_substrate = max(self.substrate.positions[:, 2])
        z_min_layer = min(self.layer.positions[:, 2])
        z_offset = z_max_substrate - z_min_layer + interface["distance"]

        return z_offset


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

    # Convert the poscar data to ASE atoms.
    substrate = poscar_to_atoms(substrate_data)
    layer = poscar_to_atoms(layer_data)

    interface_creator = InterfaceCreator(substrate, layer, SETTINGS)
    interface_structure = interface_creator.create()

    print("Interface structure: ", interface_structure)
    print("Strain along lattice a:", interface_creator.calculate_strain()["a"])
    print("Strain along lattice b:", interface_creator.calculate_strain()["b"])

    # Make the data available to the platform JS environment.
    globals()["materials_out"] = [
        {
            "poscar": atoms_to_poscar(interface_structure),
        }
    ]

    return globals()


main()
