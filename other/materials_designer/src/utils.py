import io
from pymatgen.core import Structure, Lattice
from ase import Atoms as ase_Atoms
from ase.io import read, write


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
        "constraints": [],  # Assuming there are no constraints
    }

    # Extract lattice information
    lattice = {
        "a": structure.lattice.a,
        "b": structure.lattice.b,
        "c": structure.lattice.c,
        "alpha": structure.lattice.alpha,
        "beta": structure.lattice.beta,
        "gamma": structure.lattice.gamma,
        "units": {"length": "angstrom", "angle": "degree"},
        "type": "FCC",  # You need a way to determine the lattice type
        "vectors": {
            "a": structure.lattice.matrix[0].tolist(),
            "b": structure.lattice.matrix[1].tolist(),
            "c": structure.lattice.matrix[2].tolist(),
            "alat": 1,  # This seems to be a scaling factor; adjust if necessary
            "units": "angstrom",
        },
    }

    # Combine into a material dictionary
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
    Converts pymatgen Structure object to ase.Atoms object

    Args:
        structure (Structure): pymatgen Structure object

    Returns:
        ase.Atoms: ase.Atoms object
    """
    poscar = structure.to(fmt="poscar")
    atoms = poscar_to_ase(poscar)

    return atoms
