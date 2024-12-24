import io
import numpy as np
from pymatgen.core import Structure, Lattice
from pymatgen.io.ase import AseAtomsAdaptor
from ase import Atoms
from ase.io import read, write
from ase.calculators.calculator import Calculator


# TODO: Remove this file and use Made.Tools instead
def to_pymatgen(material_data):
    """
    Converts material object in ESSE format to a pymatgen Structure object.

    Args:
        material_data (dict): A dictionary containing the material information in ESSE format.

    Returns:
        Structure: A pymatgen Structure object.
    """

    # Extract lattice information
    lattice_params = material_data["lattice"]
    a = lattice_params["a"]
    b = lattice_params["b"]
    c = lattice_params["c"]
    alpha = lattice_params["alpha"]
    beta = lattice_params["beta"]
    gamma = lattice_params["gamma"]

    # Create a Lattice from parameters
    lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)

    # Extract the basis information
    basis = material_data["basis"]
    elements = [element["value"] for element in basis["elements"]]
    coordinates = [coord["value"] for coord in basis["coordinates"]]

    # Assuming that the basis units are fractional since it's a crystal basis
    coords_are_cartesian = "units" in basis and basis["units"].lower() == "angstrom"

    # Create the Structure
    structure = Structure(lattice, elements, coordinates, coords_are_cartesian=coords_are_cartesian)

    return structure


def from_pymatgen(structure: Structure):
    """
    Converts a pymatgen Structure object to a material object in ESSE format.

    Args:
        structure (Structure): A pymatgen Structure object.

    Returns:
        dict: A dictionary containing the material information in ESSE format.
    """
    basis = {
        "elements": [{"id": i, "value": str(site.specie)} for i, site in enumerate(structure.sites)],
        "coordinates": [{"id": i, "value": list(site.frac_coords)} for i, site in enumerate(structure.sites)],
        "units": "crystal",
        "cell": structure.lattice.matrix.tolist(),
        "constraints": [],
    }

    lattice = {
        "a": structure.lattice.a,
        "b": structure.lattice.b,
        "c": structure.lattice.c,
        "alpha": structure.lattice.alpha,
        "beta": structure.lattice.beta,
        "gamma": structure.lattice.gamma,
        "units": {"length": "angstrom", "angle": "degree"},
        "type": "TRI",
        "vectors": {
            "a": structure.lattice.matrix[0].tolist(),
            "b": structure.lattice.matrix[1].tolist(),
            "c": structure.lattice.matrix[2].tolist(),
            "alat": 1,
            "units": "angstrom",
        },
    }

    material = {
        "name": structure.formula,
        "basis": basis,
        "lattice": lattice,
        "isNonPeriodic": not structure.is_ordered,
        "_id": "",
        "metadata": {"boundaryConditions": {"type": "pbc", "offset": 0}},
        "isUpdated": True,
    }

    return material


def poscar_to_ase(poscar: str):
    """
    Converts POSCAR to ase.Atoms object

    Args:
        poscar (str): POSCAR content

    Returns:
        ase.Atoms: ase.Atoms object
    """
    input_ = io.StringIO(poscar)
    atoms = read(input_, format="vasp")

    return atoms


def ase_to_poscar(atoms: Atoms):
    """
    Converts ase.Atoms object to POSCAR format

    Args:
        atoms (ase.Atoms): ase.Atoms object

    Returns:
        str: POSCAR string
    """
    output = io.StringIO()
    write(output, atoms, format="vasp")
    content = output.getvalue()
    output.close()

    return content


def ase_to_pymatgen(atoms: Atoms):
    """
    Converts ase.Atoms object to pymatgen Structure object

    Args:
        atoms (ase.Atoms): ase.Atoms object

    Returns:
        Structure: pymatgen Structure object
    """

    adaptor = AseAtomsAdaptor()
    structure = adaptor.get_structure(atoms)

    return structure


def pymatgen_to_ase(structure: Structure):
    """
    Converts a Pymatgen Structure object to an ASE Atoms object,
    transferring the 'interface_label' property to the ASE atom tags.

    Args:
        structure (Structure): The Pymatgen Structure object to convert.

    Returns:
        Atoms: The resulting ASE Atoms object with 'interface_label' as tags.
    """
    # Define a mapping from string labels to integer tags
    label_to_tag = {
        "substrate": 0,
        "film": 1,
    }

    adaptor = AseAtomsAdaptor()
    atoms = adaptor.get_atoms(structure)

    # Transfer 'interface_label' property to atom tags using the mapping
    for i, site in enumerate(structure):
        label = site.properties.get("interface_label", None)
        tag = label_to_tag.get(label, 0)  # Default to 0 if 'interface_label' is not set or is unrecognized
        atoms[i].tag = tag

    return atoms


def calculate_average_interlayer_distance(atoms, tag_substrate, tag_film, threshold=0.5):
    """
    Calculate the average distance between the top layer of substrate atoms and the bottom layer of film atoms.

    Args:
        atoms (ase.Atoms): The ASE Atoms object containing both sets of atoms.
        tag_substrate (int): The tag representing the substrate atoms.
        tag_film (int): The tag representing the film atoms.
        threshold (float): The threshold for identifying the top and bottom layers of atoms.

    Returns:
        float: The average distance between the top layer of substrate and the bottom layer of film.
    """
    # Extract z-coordinates of substrate and film atoms
    z_substrate = atoms.positions[atoms.get_tags() == tag_substrate][:, 2]
    z_film = atoms.positions[atoms.get_tags() == tag_film][:, 2]

    # Identify the top layer of substrate atoms and bottom layer of film atoms by z-coordinate
    top_substrate_layer_z = np.max(z_substrate)
    bottom_film_layer_z = np.min(z_film)

    # Get the average z-coordinate of the top substrate atoms (within a threshold from the top)
    top_substrate_atoms = z_substrate[z_substrate >= top_substrate_layer_z - threshold]
    avg_z_top_substrate = np.mean(top_substrate_atoms)

    # Get the average z-coordinate of the bottom film atoms (within a threshold from the bottom)
    bottom_film_atoms = z_film[z_film <= bottom_film_layer_z + threshold]
    avg_z_bottom_film = np.mean(bottom_film_atoms)

    # Calculate the average distance between the top layer of substrate and the bottom layer of film
    average_interlayer_distance = avg_z_bottom_film - avg_z_top_substrate
    return average_interlayer_distance


def filter_atoms_by_tag(atoms: Atoms, material_index):
    """
    Filter atoms by their tag, corresponding to the material index.

    Args:
        atoms (ase.Atoms): The Atoms object to filter.
        material_index (int): The numeric tag to filter by.

    Returns:
        ase.Atoms: The filtered Atoms object.
    """
    return atoms[atoms.get_tags() == material_index]


def get_surface_area(atoms: Atoms):
    """
    Calculate the surface area of the atoms.

    Args:
        atoms (ase.Atoms): The Atoms object to calculate the surface area of.

    Returns:
        float: The surface area of the atoms.
    """
    matrix = atoms.cell
    return np.linalg.norm(np.cross(matrix[0], matrix[1]))


def get_total_energy(atoms: Atoms, calculator: Calculator):
    """
    Set calculator for ASE Atoms and calculate the total energy.

    Args:
        atoms (ase.Atoms): The Atoms object to calculate the energy of.
        calculator (ase.calculators.calculator.Calculator): The calculator to use for the energy calculation.

    Returns:
        float: The total energy of the atoms.
    """
    atoms.set_calculator(calculator)
    return atoms.get_total_energy()


def get_total_energy_per_atom(atoms: Atoms, calculator: Calculator):
    """
        Calculate the energy per atom.

    Args:
            atoms (ase.Atoms): The Atoms object to calculate the energy of.
            calculator (ase.calculators.calculator.Calculator): The calculator to use for the energy calculation.

    Returns:
            float: The energy per atom of the atoms.
    """
    return get_total_energy(atoms, calculator) / atoms.get_global_number_of_atoms()


def get_surface_energy(slab: Atoms, bulk: Atoms, calculator: Calculator):
    """
    Calculate the surface energy by subtracting the weighted bulk energy from the slab energy.

    Args:
        slab (ase.Atoms): The slab Atoms object to calculate the surface energy of.
        bulk (ase.Atoms): The bulk Atoms object to calculate the surface energy of.
        calculator (ase.calculators.calculator.Calculator): The calculator to use for the energy calculation.

    Returns:
        float: The surface energy of the slab.
    """
    number_of_atoms = slab.get_global_number_of_atoms()
    area = get_surface_area(slab)
    return (get_total_energy(slab, calculator) - get_total_energy_per_atom(bulk, calculator) * number_of_atoms) / (
        2 * area
    )


def get_adhesion_energy(interface: Atoms, substrate_slab: Atoms, layer_slab: Atoms, calculator: Calculator):
    """
    Calculate the adhesion energy.
    The adhesion energy is the difference between the energy of the interface and the sum of the energies of the substrate and layer.
    According to: 10.1088/0953-8984/27/30/305004

    Args:
        interface (ase.Atoms): The interface Atoms object to calculate the adhesion energy of.
        substrate_slab (ase.Atoms): The substrate slab Atoms object to calculate the adhesion energy of.
        layer_slab (ase.Atoms): The layer slab Atoms object to calculate the adhesion energy of.
        calculator (ase.calculators.calculator.Calculator): The calculator to use for the energy calculation.

    Returns:
        float: The adhesion energy of the interface.
    """
    energy_substrate_slab = get_total_energy(substrate_slab, calculator)
    energy_layer_slab = get_total_energy(layer_slab, calculator)
    energy_interface = get_total_energy(interface, calculator)
    area = get_surface_area(interface)
    return (energy_substrate_slab + energy_layer_slab - energy_interface) / area


def get_interfacial_energy(
    interface: Atoms,
    substrate_slab: Atoms,
    substrate_bulk: Atoms,
    layer_slab: Atoms,
    layer_bulk: Atoms,
    calculator: Calculator,
):
    """
    Calculate the interfacial energy.
    The interfacial energy is the sum of the surface energies of the substrate and layer minus the adhesion energy.
    According to Dupré's formula

    Args:
        interface (ase.Atoms): The interface Atoms object to calculate the interfacial energy of.
        substrate_slab (ase.Atoms): The substrate slab Atoms object to calculate the interfacial energy of.
        substrate_bulk (ase.Atoms): The substrate bulk Atoms object to calculate the interfacial energy of.
        layer_slab (ase.Atoms): The layer slab Atoms object to calculate the interfacial energy of.
        layer_bulk (ase.Atoms): The layer bulk Atoms object to calculate the interfacial energy of.
        calculator (ase.calculators.calculator.Calculator): The calculator to use for the energy calculation.

    Returns:
        float: The interfacial energy of the interface.
    """

    surface_energy_substrate = get_surface_energy(substrate_slab, substrate_bulk, calculator)
    surface_energy_layer = get_surface_energy(layer_slab, layer_bulk, calculator)
    adhesion_energy = get_adhesion_energy(interface, substrate_slab, layer_slab, calculator)
    return surface_energy_layer + surface_energy_substrate - adhesion_energy


def from_ase(atoms: Atoms):
    """
    Convert an ASE Atoms object to a material object in ESSE format.

    Args:
        atoms (ase.Atoms): The ASE Atoms object to convert.

    Returns:
        dict: A dictionary containing the material information in ESSE format.
    """
    structure = ase_to_pymatgen(atoms)
    material = from_pymatgen(structure)
    return material


def translate_to_bottom(structure):
    """
    Translate the structure to the bottom of the cell.
    Args:
        structure (Structure): The pymatgen Structure object to translate.

    Returns:
        Structure: The translated pymatgen Structure object.
    """
    min_c = min(site.c for site in structure)
    translation_vector = [0, 0, -min_c]
    translated_structure = structure.copy()
    for site in translated_structure:
        site.coords += translation_vector
    return translated_structure
