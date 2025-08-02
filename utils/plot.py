from typing import Dict, List, Union

from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.interface import ZSLMatchHolder
from mat3ra.made.tools.analyze.rdf import RadialDistributionFunction
from mat3ra.utils.jupyterlite.plot import plot_distribution_function, scatter_plot_2d


def plot_strain_vs_area(matches: List["ZSLMatchHolder"], settings: Dict[str, Union[str, int]]) -> None:
    """
    Plot strain vs number of atoms for interfaces.

    Args:
        matches: List of interface matches to plot.
        settings: Plot settings.

    """
    x_values = []
    y_values = []
    hover_texts = []
    trace_names = []

    for index, match in enumerate(matches):
        strain_percentage = match.total_strain_percentage
        match_area = match.match_area

        x_values.append(strain_percentage)
        y_values.append(match_area)
        hover_texts.append(f"Index: {index}<br>Strain: {strain_percentage:.2f}%<br>Area: {match_area:.2f} Å²")
        trace_names.append(f"Index: {index}")

    plot_settings = {
        "x_title": "Strain (%)",
        "y_title": "Area (Å²)",
        "x_scale": settings["X_SCALE"],
        "y_scale": settings["Y_SCALE"],
        "height": settings["HEIGHT"],
        "legend_title": "Interfaces Indices",
    }

    fig = scatter_plot_2d(x_values, y_values, hover_texts, plot_settings, trace_names)
    fig.show()


def plot_twisted_interface_solutions(interfaces: List["Material"]) -> None:
    """
    Plot twisted interface solutions.

    Args:
        interfaces: List of interfaces to plot.
    """
    x_values = []
    y_values = []
    hover_texts = []
    trace_names = []

    for i, interface in enumerate(interfaces):
        angle = interface.metadata.get("actual_twist_angle", 0)
        size = len(interface.basis.elements.ids)

        x_values.append(angle)
        y_values.append(size)
        hover_texts.append(f"Interface {i+1}<br>Angle: {angle:.2f}°<br>Atoms: {size}<br>")
        trace_names.append(f"Interface {i+1}")

    plot_settings = {"x_title": "Twist Angle (°)", "y_title": "Number of Atoms", "title": "Twisted Interface Solutions"}

    fig = scatter_plot_2d(x_values, y_values, hover_texts, plot_settings, trace_names)
    fig.show()


def plot_rdf(material: "Material", cutoff: float = 10.0, bin_size: float = 0.1) -> None:
    """
    Plot RDF for a material.
    """
    rdf = RadialDistributionFunction.from_material(material, cutoff=cutoff, bin_size=bin_size)
    plot_distribution_function(
        rdf.bin_centers, rdf.rdf, xlabel="Distance (Å)", ylabel="g(r)", title="Radial Distribution Function (RDF)"
    )
