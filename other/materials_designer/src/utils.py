from pymatgen.core import Structure, Lattice


def to_pymatgen(material_data):
    """
    Convert material object in ESSE format to a pymatgen Structure object.

    Args:
        material_data (dict): A dictionary containing the material information in ESSE format.

    Returns:
        Structure: A pymatgen Structure object.
    """

    # Extract lattice information
    lattice_params = material_data["lattice"]
    lattice_vectors = lattice_params["vectors"]
    a = lattice_vectors["a"]
    b = lattice_vectors["b"]
    c = lattice_vectors["c"]

    # Create a Lattice
    lattice = Lattice([a, b, c])

    # Extract the basis information
    basis = material_data["basis"]
    elements = [element["value"] for element in basis["elements"]]
    coordinates = [coord["value"] for coord in basis["coordinates"]]

    # Assuming that the basis units are fractional since it's a crystal basis
    coords_are_cartesian = basis["units"].lower() != "crystal"

    # Create the Structure
    structure = Structure(lattice, elements, coordinates, coords_are_cartesian=coords_are_cartesian)

    return structure


def from_pymatgen(structure: Structure):
    """
    Convert a pymatgen Structure object to a material object in ESSE format.

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
