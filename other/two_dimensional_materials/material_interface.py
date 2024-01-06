from ase.build import surface, supercells
from ase.io import read, write
from ase.visualize import view
import io
import numpy as np


class MaterialInterface:
    def __init__(self, substrate, material, settings=None):
        self.substrate = substrate
        self.material = material
        self.settings = globals()["default_settings"]
        if settings:
            for key in self.settings.keys():
                if key in settings:
                    self.settings[key].update(settings[key])
        self.structure = self.create_structure()

    def create_structure(self):
        self.materials = []
        slab = self.settings["slab"]
        interface = self.settings["interface"]

        self.substrate = surface(
            self.substrate,
            (slab["miller:h"], slab["miller:k"], slab["miller:l"]),
            vacuum=slab["vacuum"],
            layers=slab["number_of_layers"],
        )

        slab_v_matrix = self.expand_matrix_2x2_to_3x3(interface["slab_v:matrix"])
        layer_v_matrix = self.expand_matrix_2x2_to_3x3(interface["layer_v:matrix"])

        self.substrate = supercells.make_supercell(self.substrate, slab_v_matrix)
        # 0
        self.materials.append(self.substrate.copy())

        vector1 = self.substrate.get_cell()[0]
        vector2 = self.material.get_cell()[0]

        angle = get_angle(vector1, vector2)
        print(angle)
        self.material.rotate(-angle, "z", rotate_cell=True)
        # 1
        self.materials.append(self.material.copy())

        self.substrate.wrap()
        # 2
        self.materials.append(self.substrate.copy())

        self.material = supercells.make_supercell(self.material, layer_v_matrix)
        # 3
        self.materials.append(self.material.copy())

        self.original_material = self.material.copy()

        self.material.set_cell(self.substrate.get_cell(), scale_atoms=True)
        # 4
        self.materials.append(self.material.copy())

        self.material.wrap()
        # 5
        self.materials.append(self.material.copy())

        z_max_substrate = max(self.substrate.positions[:, 2])
        z_min_material = min(self.material.positions[:, 2])
        z_offset = z_max_substrate - z_min_material + interface["distance"]
        self.material.positions[:, 2] += z_offset

        interface = self.substrate + self.material
        # 6
        self.materials.append(interface.copy())
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

    def view(self, material=None, repeat=(3, 3, 1)):
        if material is None:
            material = self.structure
        view(material * repeat)

    def poscar_to_atoms(self, poscar):
        input = io.StringIO(poscar)
        atoms = read(input, format="vasp")
        return atoms

    def write_atoms_to_poscar(atoms):
        output = io.StringIO()
        write(output, atoms, format="vasp")
        content = output.getvalue()
        output.close()
        return content

    def poscar_to_cif(poscar):
        input = io.StringIO(poscar)
        atoms = read(input, format="vasp")

        output = io.BytesIO()
        write(output, atoms, format="cif")
        content = output.getvalue()
        output.close()
        return content

    def get_niggli_cell(poscar):
        structure = Structure.from_str(poscar, "poscar")
        structure = mg.symmetry.analyzer.SpacegroupAnalyzer(structure).get_primitive_standard_structure()
        struct_poscar = structure.to(fmt="poscar")
        return struct_poscar

    def expand_matrix_2x2_to_3x3(matrix_2x2):
        matrix_3x3 = [[0, 0, 0], [0, 0, 0], [0, 0, 1]]

        for i in range(2):
            for j in range(2):
                matrix_3x3[i][j] = matrix_2x2[i][j]

        return matrix_3x3

    def dotproduct(v1, v2):
        return sum((a * b) for a, b in zip(v1, v2))

    def length(v):
        return math.sqrt(dotproduct(v, v))

    def get_angle(v1, v2):
        return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2))) * 180 / math.pi
