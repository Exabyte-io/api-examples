"""TITLE: Strain Matching of a 2D Material Layer on a Surface."""

"""BLOCK: Install Dependencies"""

"""
This scripts installs the required packages for a basic Pyodide environment for Materials structure manipulation.
It uses micropip - a Python package manager for Pyodide.
"""
import micropip

packages = [
    "https://files.mat3ra.com:44318/uploads/pymatgen-2023.9.10-py3-none-any.whl",
    "https://files.mat3ra.com:44318/web/pyodide/spglib-2.0.2-py3-none-any.whl",
    "https://files.pythonhosted.org/packages/d9/0e/2a05efa11ea33513fbdf4a2e2576fe94fd8fa5ad226dbb9c660886390974/ruamel.yaml-0.17.32-py3-none-any.whl",
    "ase==3.22.1",
    "networkx==3.2.1",
    "monty==2023.11.3",
    "scipy==1.11.1",
    "lzma",
    "tabulate==0.9.0",
    "sqlite3",
    "sympy==1.12",
    "uncertainties==3.1.6",
]


async def install_package(pkg):
    """
    Installs a package in a Pyodide environment.
    Args:
        pkg: The name of the package to install.

    Returns:
        None
    """
    is_url = pkg.startswith("http://") or pkg.startswith("https://")
    are_dependencies_installed = not is_url
    await micropip.install(pkg, deps=are_dependencies_installed)
    # Extract package name for printing
    pkg_name = pkg.split("/")[-1].split("-")[0] if is_url else pkg.split("==")[0]
    print(f"Installed {pkg_name}")


for package in packages:
    await install_package(package)

"""BLOCK: Utils, Class Definitions, and main()"""
"""
This script employs the Zur and McGill SuperLattice algorithm (Journal of Applied Physics 55 (1984), 378 ; doi: 10.1063/1.333084)  for strain matching a 2D material layer on a surface. 
Using pymatgen, it constructs coherent interfaces between a substrate and a film layer, with a focus on lattice matching and interface terminations. 
Link to pymatgen ZSL implementation: https://pymatgen.org/pymatgen.analysis.interfaces.html#pymatgen.analysis.interfaces.zsl.ZSLGenerator
Key parameters like Miller indices and layer thicknesses can be customized.
The result of the lagorithm is a list of interfaces, sorted by the mean absolute strain.
Plot shows the mean absolute strain vs. the number of atoms in the interface to select corresponding material by index. 
"""
from __future__ import annotations


# Select materials from the list of input materials under `materials_in` in globals().
SUBSTRATE_INDEX = 0
LAYER_INDEX = 1

# Select interfaces from list of generated interfaces sorted by increasing mean absolute strain
# By interval or by index, e.g. [0, 10] or 0
OUTPUT_INDICES = [0, 10]

# Select Miller indices and thickness for substrate and layer.
SUBSTRATE_MILLER = (1, 1, 1)
SUBSTRATE_THICKNESS = 3
LAYER_MILLER = (0, 0, 1)
LAYER_THICKNESS = 1

# Select distance between layers
DISTANCE = 3.0

# Select parameters for the ZSL algorithm
# as defined in: https://pymatgen.org/pymatgen.analysis.interfaces.html#pymatgen.analysis.interfaces.zsl.ZSLGenerator
MAX_AREA = 400
MAX_AREA_TOL = 0.09
MAX_LENGTH_TOL = 0.03
MAX_ANGLE_TOL = 0.01

# Strains within this tolerance are considered equal
STRAIN_TOL = 10e-6

# Plot settings
# Strain axis limits in percent
X_MIN = 0.01
X_MAX = 100

# Number of atoms axis limits
Y_MIN = 1
Y_MAX = 1000


"""
NOTE: DO NOT edit code below unless you know what you are doing.
"""

""" 
Classes and Definitions
"""

from itertools import product
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
from numpy.testing import assert_allclose
from scipy.linalg import polar

import pymatgen
from pymatgen.analysis.interfaces.zsl import ZSLGenerator as PymatgenZSL
from pymatgen.core.structure import Structure
from operator import itemgetter
from pymatgen.analysis.elasticity.strain import Deformation
from pymatgen.analysis.interfaces.zsl import ZSLGenerator, fast_norm
from pymatgen.core.interface import Interface, label_termination
from pymatgen.core.surface import SlabGenerator

if TYPE_CHECKING:
    from collections.abc import Iterator, Sequence

    from pymatgen.core import Structure


class CoherentInterfaceBuilder:
    """
    This class constructs the coherent interfaces between two crystalline slabs
    Coherency is defined by matching lattices not sub-planes.
    """

    def __init__(
        self,
        substrate_structure: Structure,
        film_structure: Structure,
        film_miller: tuple[int, int, int],
        substrate_miller: tuple[int, int, int],
        zslgen: ZSLGenerator | None = None,
    ):
        """
        Args:
            substrate_structure: structure of substrate
            film_structure: structure of film
            film_miller: miller index of the film layer
            substrate_miller: miller index for the substrate layer
            zslgen: BiDirectionalZSL if you want custom lattice matching tolerances for coherency.
        """
        # Bulk structures
        self.substrate_structure = substrate_structure
        self.film_structure = film_structure
        self.film_miller = film_miller
        self.substrate_miller = substrate_miller
        self.zslgen = zslgen or ZSLGenerator(bidirectional=True)

        self._find_matches()
        self._find_terminations()

    def _find_matches(self) -> None:
        """
        Finds and stores the ZSL matches.
        """
        self.zsl_matches = []

        film_sg = SlabGenerator(
            self.film_structure,
            self.film_miller,
            min_slab_size=1,
            min_vacuum_size=3,
            in_unit_planes=True,
            center_slab=True,
            primitive=True,
            reorient_lattice=False,  # This is necessary to not screw up the lattice
        )

        sub_sg = SlabGenerator(
            self.substrate_structure,
            self.substrate_miller,
            min_slab_size=1,
            min_vacuum_size=3,
            in_unit_planes=True,
            center_slab=True,
            primitive=True,
            reorient_lattice=False,  # This is necessary to not screw up the lattice
        )

        film_slab = film_sg.get_slab(shift=0)
        sub_slab = sub_sg.get_slab(shift=0)

        film_vectors = film_slab.lattice.matrix
        substrate_vectors = sub_slab.lattice.matrix

        # Generate all possible interface matches
        self.zsl_matches = list(self.zslgen(film_vectors[:2], substrate_vectors[:2], lowest=False))

        for match in self.zsl_matches:
            xform = get_2d_transform(film_vectors, match.film_vectors)
            strain, rot = polar(xform)
            assert_allclose(
                strain, np.round(strain), atol=1e-12
            ), "Film lattice vectors changed during ZSL match, check your ZSL Generator parameters"

            xform = get_2d_transform(substrate_vectors, match.substrate_vectors)
            strain, rot = polar(xform)
            assert_allclose(
                strain, strain.astype(int), atol=1e-12
            ), "Substrate lattice vectors changed during ZSL match, check your ZSL Generator parameters"

    def _find_terminations(self):
        """
        Finds all terminations.
        """
        film_sg = SlabGenerator(
            self.film_structure,
            self.film_miller,
            min_slab_size=1,
            min_vacuum_size=3,
            in_unit_planes=True,
            center_slab=True,
            primitive=True,
            reorient_lattice=False,  # This is necessary to not screw up the lattice
        )

        sub_sg = SlabGenerator(
            self.substrate_structure,
            self.substrate_miller,
            min_slab_size=1,
            min_vacuum_size=3,
            in_unit_planes=True,
            center_slab=True,
            primitive=True,
            reorient_lattice=False,  # This is necessary to not screw up the lattice
        )

        film_slabs = film_sg.get_slabs()
        sub_slabs = sub_sg.get_slabs()

        film_shits = [s.shift for s in film_slabs]
        film_terminations = [label_termination(s) for s in film_slabs]

        sub_shifts = [s.shift for s in sub_slabs]
        sub_terminations = [label_termination(s) for s in sub_slabs]

        self._terminations = {
            (film_label, sub_label): (film_shift, sub_shift)
            for (film_label, film_shift), (sub_label, sub_shift) in product(
                zip(film_terminations, film_shits), zip(sub_terminations, sub_shifts)
            )
        }
        self.terminations = list(self._terminations)

    def get_interfaces(
        self,
        termination: tuple[str, str],
        gap: float = 2.0,
        vacuum_over_film: float = 20.0,
        film_thickness: float = 1,
        substrate_thickness: float = 1,
        in_layers: bool = True,
    ) -> Iterator[Interface]:
        """
        Generates interface structures given the film and substrate structure
        as well as the desired terminations.

        Args:
            termination (tuple[str, str]): termination from self.termination list
            gap (float, optional): gap between film and substrate. Defaults to 2.0.
            vacuum_over_film (float, optional): vacuum over the top of the film. Defaults to 20.0.
            film_thickness (float, optional): the film thickness. Defaults to 1.
            substrate_thickness (float, optional): substrate thickness. Defaults to 1.
            in_layers (bool, optional): set the thickness in layer units. Defaults to True.

        Yields:
            Iterator[Interface]: interfaces from slabs
        """
        film_sg = SlabGenerator(
            self.film_structure,
            self.film_miller,
            min_slab_size=film_thickness,
            min_vacuum_size=3,
            in_unit_planes=in_layers,
            center_slab=True,
            primitive=True,
            reorient_lattice=False,  # This is necessary to not screw up the lattice
        )

        sub_sg = SlabGenerator(
            self.substrate_structure,
            self.substrate_miller,
            min_slab_size=substrate_thickness,
            min_vacuum_size=3,
            in_unit_planes=in_layers,
            center_slab=True,
            primitive=True,
            reorient_lattice=False,  # This is necessary to not screw up the lattice
        )

        film_shift, sub_shift = self._terminations[termination]

        film_slab = film_sg.get_slab(shift=film_shift)
        sub_slab = sub_sg.get_slab(shift=sub_shift)

        for match in self.zsl_matches:
            # Build film superlattice
            super_film_transform = np.round(
                from_2d_to_3d(get_2d_transform(film_slab.lattice.matrix[:2], match.film_sl_vectors))
            ).astype(int)
            film_sl_slab = film_slab.copy()
            film_sl_slab.make_supercell(super_film_transform)
            assert_allclose(
                film_sl_slab.lattice.matrix[2], film_slab.lattice.matrix[2], atol=1e-08
            ), "2D transformation affected C-axis for Film transformation"
            assert_allclose(
                film_sl_slab.lattice.matrix[:2], match.film_sl_vectors, atol=1e-08
            ), "Transformation didn't make proper supercell for film"

            # Build substrate superlattice
            super_sub_transform = np.round(
                from_2d_to_3d(get_2d_transform(sub_slab.lattice.matrix[:2], match.substrate_sl_vectors))
            ).astype(int)
            sub_sl_slab = sub_slab.copy()
            sub_sl_slab.make_supercell(super_sub_transform)
            assert_allclose(
                sub_sl_slab.lattice.matrix[2], sub_slab.lattice.matrix[2], atol=1e-08
            ), "2D transformation affected C-axis for Film transformation"
            assert_allclose(
                sub_sl_slab.lattice.matrix[:2], match.substrate_sl_vectors, atol=1e-08
            ), "Transformation didn't make proper supercell for substrate"

            # Add extra info
            match_dict = match.as_dict()
            interface_properties = {k: match_dict[k] for k in match_dict if not k.startswith("@")}

            dfm = Deformation(match.match_transformation)

            strain = dfm.green_lagrange_strain
            interface_properties["strain"] = strain
            interface_properties["von_mises_strain"] = strain.von_mises_strain
            interface_properties["termination"] = termination
            interface_properties["film_thickness"] = film_thickness
            interface_properties["substrate_thickness"] = substrate_thickness

            yield {
                "interface": Interface.from_slabs(
                    substrate_slab=sub_sl_slab,
                    film_slab=film_sl_slab,
                    gap=gap,
                    vacuum_over_film=vacuum_over_film,
                    interface_properties=interface_properties,
                ),
                "strain": strain,
                "von_mises_strain": strain.von_mises_strain,
                "mean_abs_strain": round(np.mean(np.abs(strain)) / STRAIN_TOL) * STRAIN_TOL,
                "film_sl_vectors": match.film_sl_vectors,
                "substrate_sl_vectors": match.substrate_sl_vectors,
                "film_transform": super_film_transform,
                "substrate_transform": super_sub_transform,
            }


def get_rot_3d_for_2d(film_matrix, sub_matrix) -> np.ndarray:
    """
    Find transformation matrix that will rotate and strain the film to the substrate while preserving the c-axis.
    """
    film_matrix = np.array(film_matrix)
    film_matrix = film_matrix.tolist()[:2]
    film_matrix.append(np.cross(film_matrix[0], film_matrix[1]))

    # Generate 3D lattice vectors for substrate super lattice
    # Out of plane substrate super lattice has to be same length as
    # Film out of plane vector to ensure no extra deformation in that
    # direction
    sub_matrix = np.array(sub_matrix)
    sub_matrix = sub_matrix.tolist()[:2]
    temp_sub = np.cross(sub_matrix[0], sub_matrix[1]).astype(float)  # conversion to float necessary if using numba
    temp_sub = temp_sub * fast_norm(
        np.array(film_matrix[2], dtype=float)
    )  # conversion to float necessary if using numba
    sub_matrix.append(temp_sub)

    transform_matrix = np.transpose(np.linalg.solve(film_matrix, sub_matrix))

    rot, _ = polar(transform_matrix)

    return rot


def get_2d_transform(start: Sequence, end: Sequence) -> np.ndarray:
    """
    Gets a 2d transformation matrix
    that converts start to end.
    """
    return np.dot(end, np.linalg.pinv(start))


def from_2d_to_3d(mat: np.ndarray) -> np.ndarray:
    """
    Converts a 2D matrix to a 3D matrix.
    """
    new_mat = np.diag([1.0, 1.0, 1.0])
    new_mat[:2, :2] = mat
    return new_mat


def plot_strain_vs_atoms(strain_mode, sorted_interfaces):
    """
    Plots the strain vs. the number of atoms in the interface. With hover-over labels.
    """
    fig, ax = plt.subplots()

    # Scatter plot
    x = [i[strain_mode] * 100 for i in sorted_interfaces]  # in precentage
    y = [i["interface"].num_sites for i in sorted_interfaces]
    sc = ax.scatter(x, y)

    # Annotation for the hover-over labels
    annot = ax.annotate(
        "",
        xy=(0, 0),
        xytext=(20, 20),
        textcoords="offset points",
        bbox=dict(boxstyle="round", fc="w"),
        arrowprops=dict(arrowstyle="->"),
    )
    annot.set_visible(False)

    def update_annot(ind):
        pos = sc.get_offsets()[ind["ind"][0]]
        annot.xy = pos
        text = "{}".format(" ".join([str(index) for index in ind["ind"]]))
        annot.set_text(text)
        annot.get_bbox_patch().set_alpha(0.4)

    def hover(event):
        vis = annot.get_visible()
        if event.inaxes == ax:
            cont, ind = sc.contains(event)
            if cont:
                update_annot(ind)
                annot.set_visible(True)
                fig.canvas.draw_idle()
            else:
                if vis:
                    annot.set_visible(False)
                    fig.canvas.draw_idle()

    # Connect the hover event
    fig.canvas.mpl_connect("motion_notify_event", hover)

    # Set the scale and labels
    plt.xscale("log")
    plt.yscale("log")
    plt.xlim(0.01, 100)
    plt.ylim(1, 1000)
    plt.xlabel(strain_mode)
    plt.ylabel("number of atoms")

    plt.show()


def output_materials(sorted_interfaces, output_indices, strain_mode):
    """
    Outputs the materials in the output range.
    """
    # Initialize an empty list to store the materials
    materials_out = []

    # Check if output_indices is a single number or a list
    if isinstance(output_indices, int):
        # If it's a single number, output the material with that index
        output_range = [output_indices]
    else:
        # If it's a list, output materials in that range
        output_range = range(*output_indices)

    # Loop over the interfaces in the output range
    for index in output_range:
        best_interface = sorted_interfaces[index]["interface"]

        wrapped_structure = best_interface.copy()
        # This should wrap the atoms inside the unit cell
        wrapped_structure.make_supercell([1, 1, 1])
        interface_poscar = wrapped_structure.to(fmt="poscar")

        # Add the strain to the POSCAR first line to use as a name
        strain = sorted_interfaces[index][strain_mode]
        lines = interface_poscar.split("\n")
        lines[0] = lines[0] + " strain: " + "{:.3f}".format(strain * 100) + "%"
        interface_poscar = "\n".join(lines)

        # Append the material to the list
        materials_out.append({"poscar": interface_poscar})

    return materials_out


""" MAIN """


def main():
    # Interaction with Platform
    materials_in = globals()["materials_in"]
    substrate = Structure.from_str(materials_in[SUBSTRATE_INDEX].getAsPOSCAR(), fmt="poscar")
    layer = Structure.from_str(materials_in[LAYER_INDEX].getAsPOSCAR(), fmt="poscar")

    # Create Interface Builder class
    zsl = PymatgenZSL(
        max_area_ratio_tol=MAX_AREA_TOL, max_area=MAX_AREA, max_length_tol=MAX_LENGTH_TOL, max_angle_tol=MAX_ANGLE_TOL
    )
    cib = CoherentInterfaceBuilder(
        substrate_structure=substrate,
        film_structure=layer,
        substrate_miller=SUBSTRATE_MILLER,
        film_miller=LAYER_MILLER,
        zslgen=zsl,
    )

    # Run the Interface Building process
    cib._find_terminations()
    matches = cib.zsl_matches
    terminations = cib.terminations

    # Create interfaces
    interfaces = []
    for termination in terminations:
        interfaces = list(
            cib.get_interfaces(
                termination,
                gap=DISTANCE,
                film_thickness=LAYER_THICKNESS,
                substrate_thickness=SUBSTRATE_THICKNESS,
                in_layers=True,
            )
        )

    print(f"Found {len(matches)} interfaces")
    print(f"Found {len(terminations)} terminations:", terminations)

    strain_modes = {"VON_MISES": "von_mises_strain", "STRAIN": "strain", "MEAN": "mean_abs_strain"}
    strain_mode = strain_modes["MEAN"]
    interfaces_list = list(interfaces)

    # Sort interfaces by ascending strain and then by ascending number of atoms
    sorted_interfaces = sorted(interfaces_list, key=lambda x: (itemgetter(strain_mode)(x), x["interface"].num_sites))

    print("Interface with lowest strain (index 0):")
    print("    strain:", sorted_interfaces[0][strain_mode] * 100, "%")
    print("    number of atoms:", sorted_interfaces[0]["interface"].num_sites)

    # plot stran vs number of atoms via matplotlib
    plot_strain_vs_atoms(strain_mode, sorted_interfaces)

    # Return created materials to the platform
    globals()["materials_out"] = output_materials(sorted_interfaces, OUTPUT_INDICES, strain_mode)

    return globals()


main()