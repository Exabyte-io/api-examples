# data that comes from JS, here are some examples
poscars = {
    "Ni": """Ni4
1.0
   3.4751458659480110    0.0000000000000000    0.0000000000000002
   0.0000000000000006    3.4751458659480110    0.0000000000000002
   0.0000000000000000    0.0000000000000000    3.4751458659480110
Ni
4
direct
   0.0000000000000000    0.0000000000000000    0.0000000000000000 Ni
   0.0000000000000000    0.5000000000000000    0.5000000000000000 Ni
   0.5000000000000000    0.0000000000000000    0.5000000000000000 Ni
   0.5000000000000000    0.5000000000000000    0.0000000000000000 Ni
""",
    "Cu": """Cu4
1.0
   3.5774306715697510    0.0000000000000000    0.0000000000000002
   0.0000000000000006    3.5774306715697510    0.0000000000000002
   0.0000000000000000    0.0000000000000000    3.5774306715697510
Cu
4
direct
   0.0000000000000000    0.0000000000000000    0.0000000000000000 Cu
   0.0000000000000000    0.5000000000000000    0.5000000000000000 Cu
   0.5000000000000000    0.0000000000000000    0.5000000000000000 Cu
   0.5000000000000000    0.5000000000000000    0.0000000000000000 Cu
""",
    "Au": """Au4
1.0
   4.1712885314747270    0.0000000000000000    0.0000000000000003
   0.0000000000000007    4.1712885314747270    0.0000000000000003
   0.0000000000000000    0.0000000000000000    4.1712885314747270
Au
4
direct
   0.0000000000000000    0.0000000000000000    0.0000000000000000 Au
   0.0000000000000000    0.5000000000000000    0.5000000000000000 Au
   0.5000000000000000    0.0000000000000000    0.5000000000000000 Au
   0.5000000000000000    0.5000000000000000    0.0000000000000000 Au
""",
    "SiC": """Si4 C4
1.0
   4.3539932475828609    0.0000000000000000    0.0000000000000003
   0.0000000000000007    4.3539932475828609    0.0000000000000003
   0.0000000000000000    0.0000000000000000    4.3539932475828609
Si C
4 4
direct
   0.7500000000000000    0.2500000000000000    0.7500000000000000 Si4+
   0.7500000000000000    0.7500000000000000    0.2500000000000000 Si4+
   0.2500000000000000    0.2500000000000000    0.2500000000000000 Si4+
   0.2500000000000000    0.7500000000000000    0.7500000000000000 Si4+
   0.0000000000000000    0.0000000000000000    0.0000000000000000 C4-
   0.0000000000000000    0.5000000000000000    0.5000000000000000 C4-
   0.5000000000000000    0.0000000000000000    0.5000000000000000 C4-
   0.5000000000000000    0.5000000000000000    0.0000000000000000 C4-
""",
    "Graphene": """Graphene
1.0
   2.467291000	   0.000000000	   0.000000000
  -1.233645000	   2.136737000	   0.000000000
   0.000000000	   0.000000000	   7.803074000
C
2
direct
   0.000000000    0.000000000    0.000000000  C
   0.333333000    0.666667000    0.000000000  C
""",
}

# definition of MaterialInterface class and helper functions for it
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
            print("ASE is not installed, cannot view the structure")
            raise
     

# actual Pyodide code

global settings, poscar_data
poscar_data = poscars


def func():
    """This function is a gateway to Pyodide in Materials Designer"""

    poscar_data = globals()["poscar_data"]
    settings = globals()["settings"]
    materials = {}

    for material_name, poscar in poscar_data.items():
        materials[material_name] = poscar_to_atoms(poscar)

    substrate = materials["Au"]
    material = materials["Graphene"]
    interface = MaterialInterface(substrate, material, settings)
    write("structure.poscar", interface.structure, format="vasp")
    print(interface.structure)
    interface.view(repeat=(2, 2, 1))
    print("strain (a, b):", interface.calculate_strain())


func()
