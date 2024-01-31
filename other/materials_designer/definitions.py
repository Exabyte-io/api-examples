from IPython.display import display, Javascript
import json


def submit(materials):
    python_data = materials
    serialized_data = json.dumps({"materials": python_data})
    js_code = f"""
    (function() {{
        window.sendDataToHost({serialized_data})
        console.log({serialized_data})
    }})();
    """

    display(Javascript(js_code))
    print("materials sent")


def get_materials():
    js_code = """
    (function() {
        if (window.requestDataFromHost) {
            window.requestDataFromHost();
            
        } else {
            console.error('requestDataFromHost function is not defined on the window object.');
        }
    })();
    """

    display(Javascript(js_code))
    print("materials requested")


from pymatgen.core import Structure, Lattice


def to_pymatgen(material_data):
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
    # wrap atoms into the unit cell
    structure = structure.to_unit_cell()

    return structure


def from_pymatgen(structure: Structure):
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
        "metadata": {"boundaryConditions": {"type": "bc2", "offset": 0}},
        "isUpdated": True,
    }

    return material
