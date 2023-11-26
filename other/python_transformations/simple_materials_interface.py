from ase.build import surface, supercells
from ase.io import read, write
import io
import numpy as np


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


# default values
globals().setdefault("data_in", {"materials": [{"poscar": ""}, {"poscar": ""}]})

# Set the parameters
globals()["data_in"]["settings"] = {
    "slab": {
        "miller:h": 1,
        "miller:k": 1,
        "miller:l": 1,
        "vacuum": 5,
        "number_of_layers": 3,
    },
    "interface": {"slab_v:matrix": [[1, 0], [0, 1]], "layer_v:matrix": [[1, 0], [0, 1]], "distance": 2.0},
}


# Run the interface creation


def func():
    """This function is a gateway to Pyodide in Materials Designer"""
    try:
        settings = globals()["data_in"]["settings"]
        materials = globals()["data_in"]["materials"]
        substrate_data = materials[0]
        layer_data = materials[1]

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
