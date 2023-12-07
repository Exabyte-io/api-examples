"""Simple Materials Interface"""

"""BLOCK: Packages Import"""
import micropip

await micropip.install("ase")
from ase.build import surface, supercells
from ase.io import read, write
import io
import numpy as np

print("installed ase")


"""BLOCK: Classes and Definitions"""
# Parameters of the interface
SUBSTRATE_INDEX = 0
LAYER_INDEX = 1

SLAB_MILLER_H = 1
SLAB_MILLER_K = 1
SLAB_MILLER_L = 1
SLAB_VACUUM = 5
SLAB_NUMBER_OF_LAYERS = 3

INTERFACE_SLAB_V_MATRIX = [[1, 0], [0, 1]]
INTERFACE_LAYER_V_MATRIX = [[1, 0], [0, 1]]
INTERFACE_DISTANCE = 2.0


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


class MaterialInterface:
    def __init__(self, substrate, material, settings=None):
        self.substrate = substrate
        self.material = material
        self.settings = globals()["data_in"]["settings"]
        if settings:
            for key in self.settings.keys():
                if key in settings:
                    self.settings[key].update(settings[key])
        self.structure = self.create_structure()

    def create_structure(self):
        slab = self.settings["slab"]
        interface = self.settings["interface"]

        self.substrate = surface(
            self.substrate,
            (slab["miller:h"], slab["miller:k"], slab["miller:l"]),
            vacuum=slab["vacuum"],
            layers=slab["number_of_layers"],
        )

        slab_v_matrix = expand_matrix_2x2_to_3x3(interface["slab_v:matrix"])
        layer_v_matrix = expand_matrix_2x2_to_3x3(interface["layer_v:matrix"])

        self.substrate = supercells.make_supercell(self.substrate, slab_v_matrix)
        self.substrate.wrap()
        self.material = supercells.make_supercell(self.material, layer_v_matrix)
        self.original_material = self.material.copy()
        # self.material.set_cell(self.substrate.get_cell(), scale_atoms=True)
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


# Function that gets executed
def func():
    """This function gets executed and returns transformed materials to platform JS environment"""
    globals().setdefault("data_in", {"materials": [{"poscar": ""}, {"poscar": ""}]})
    globals()["data_in"]["settings"] = {
        "slab": {
            "miller:h": SLAB_MILLER_H,
            "miller:k": SLAB_MILLER_K,
            "miller:l": SLAB_MILLER_L,
            "vacuum": SLAB_VACUUM,
            "number_of_layers": SLAB_NUMBER_OF_LAYERS,
        },
        "interface": {
            "slab_v:matrix": INTERFACE_SLAB_V_MATRIX,
            "layer_v:matrix": INTERFACE_LAYER_V_MATRIX,
            "distance": INTERFACE_DISTANCE,
        },
    }
    try:
        settings = globals()["data_in"]["settings"]
        materials = globals()["data_in"]["materials"]
        substrate_data = materials[SUBSTRATE_INDEX]
        layer_data = materials[LAYER_INDEX]

        substrate = ase_poscar_to_atoms(substrate_data["poscar"])
        layer = ase_poscar_to_atoms(layer_data["poscar"])

        interface = MaterialInterface(substrate, layer, settings)

        print("Interface: ", interface.structure)
        print("strain (a, b):", interface.calculate_strain())

        globals()["data_out"]["materials"] = [
            {
                "poscar": ase_atoms_to_poscar(interface.structure),
                "metadata": {},
            }
        ]
    except Exception as e:
        print(e)

    return globals()


func()
