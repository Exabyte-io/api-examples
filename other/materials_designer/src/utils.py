import io
from pymatgen.core import Structure, Lattice
from pymatgen.io.ase import AseAtomsAdaptor
from ase import Atoms as ase_Atoms
from ase.io import read, write
import numpy as np

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
    coords_are_cartesian = 'units' in basis and basis['units'].lower() == 'angstrom'

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


def ase_to_poscar(atoms: ase_Atoms):
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


def ase_to_pymatgen(atoms: ase_Atoms):
    """
    Converts ase.Atoms object to pymatgen Structure object

    Args:
        atoms (ase.Atoms): ase.Atoms object

    Returns:
        Structure: pymatgen Structure object
    """
    poscar = ase_to_poscar(atoms)
    structure = Structure.from_str(poscar, fmt="poscar")

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
        label = site.properties.get('interface_label', None)
        tag = label_to_tag.get(label, 0)  # Default to 0 if 'interface_label' is not set or is unrecognized
        atoms[i].tag = tag

    return atoms


def calculate_average_interlayer_distance(atoms, tag_substrate, tag_film, threshold = 0.5):
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
