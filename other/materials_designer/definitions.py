from IPython.display import display, Javascript
import json


def set_materials(materials):
    """
    This function takes a Python object, serializes it to JSON, and sends it to the host environment
    through a JavaScript function defined in the JupyterLite extension `data_bridge`.

    Args:
        materials (object): The Python object to be sent to the host environment.
    """
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
    """
    This function requests materials from the host environment through a JavaScript function defined in the JupyterLite
    extension `data_bridge`. The materials are then returned to the Python environment.
    """
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
        "metadata": {"boundaryConditions": {"type": "bc2", "offset": 0}},
        "isUpdated": True,
    }

    return material


import matplotlib.pyplot as plt


def plot_strain_vs_atoms(strain_mode, sorted_interfaces, settings):
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
    plt.xscale(settings["X_SCALE"])
    plt.yscale(settings["Y_SCALE"])
    plt.xlim(settings["X_MIN"], settings["X_MAX"])
    plt.ylim(settings["Y_MIN"], settings["Y_MAX"])

    plt.xlabel("strain in %")
    plt.ylabel("number of atoms")

    plt.show()
