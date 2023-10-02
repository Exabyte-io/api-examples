from ase.build import surface, supercells
from ase.io import read, write
import io
import numpy as np

globals()["settings"] = {
    "slab": {
        "miller:h": 1,
        "miller:k": 1,
        "miller:l": 1,
        "vacuum": 1,
        "number_of_layers": 3,
    },
    "interface": {"slab_v:matrix": [[1, 0], [0, 1]], "layer_v:matrix": [[1, 0], [0, 1]], "distance": 2.0},
}


def poscar_to_atoms(poscar):
    input = io.StringIO(poscar)
    atoms = read(input, format="vasp")
    return atoms


def write_atoms_to_poscar(atoms):
    output = io.StringIO()
    write(output, atoms, format="vasp")
    content = output.getvalue()
    output.close()
    return content


def expand_matrix_2x2_to_3x3(matrix_2x2):
    matrix_3x3 = [[0, 0, 0], [0, 0, 0], [0, 0, 1]]

    for i in range(2):
        for j in range(2):
            matrix_3x3[i][j] = matrix_2x2[i][j]

    return matrix_3x3


class MaterialInterface:
    def __init__(self, substrate, material, settings=None):
        self.substrate = substrate
        self.material = material
        self.settings = globals()["settings"]
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
        self.material.set_cell(self.substrate.get_cell(), scale_atoms=True)
        self.material.wrap()

        z_offset = self.calculate_distance()
        self.material.positions[:, 2] += z_offset

        return self.substrate + self.material

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

    def view(self, material=None, repeat=(1, 1, 1)):
        if material is None:
            material = self.structure
        try :
            from ase.visualize import view    
            view(material * repeat)
        except:
            print("ASE.gui is not installed, cannot view the structure")
            raise
     

# Running the interface creation

def func():
    """This function is a gateway to Pyodide in Materials Designer"""

    poscar_data = globals()["poscar_data"]
    settings = globals()["settings"]
    substrate = poscar_to_atoms(poscar_data[0])
    material = poscar_to_atoms(poscar_data[1])

    interface = MaterialInterface(substrate, material, settings)

    print(interface.structure)
    print("strain (a, b):", interface.calculate_strain())
    return interface

func()
