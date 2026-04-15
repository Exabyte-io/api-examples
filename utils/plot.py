from typing import Dict, List, Tuple, Union

from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.interface import ZSLMatchHolder
from mat3ra.made.tools.analyze.rdf import RadialDistributionFunction
from mat3ra.utils.jupyterlite.plot import plot_distribution_function, render_figure, scatter_plot_2d
from matplotlib import pyplot as plt


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
    render_figure(fig)


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
        hover_texts.append(f"Interface {i + 1}<br>Angle: {angle:.2f}°<br>Atoms: {size}<br>")
        trace_names.append(f"Interface {i + 1}")

    plot_settings = {"x_title": "Twist Angle (°)", "y_title": "Number of Atoms", "title": "Twisted Interface Solutions"}

    fig = scatter_plot_2d(x_values, y_values, hover_texts, plot_settings, trace_names)
    render_figure(fig)


def plot_rdf(material: "Material", cutoff: float = 10.0, bin_size: float = 0.1) -> None:
    """
    Plot RDF for a material.
    """
    rdf = RadialDistributionFunction.from_material(material, cutoff=cutoff, bin_size=bin_size)
    plot_distribution_function(
        rdf.bin_centers, rdf.rdf, xlabel="Distance (Å)", ylabel="g(r)", title="Radial Distribution Function (RDF)"
    )


def plot_series(
    series: List[Dict],
    x_key: str,
    y_key: str,
    xlabel: str,
    ylabel: str,
    title: str,
    figsize: Tuple[int, int] = (8, 5),
    marker: str = "o",
    rotation: int = 45,
) -> None:
    """
    Plot a series of data points with configurable parameters.

    Args:
        series: List of dictionaries containing data points.
        x_key: Key to extract x values from series items.
        y_key: Key to extract y values from series items.
        xlabel: Label for x-axis.
        ylabel: Label for y-axis.
        title: Title of the plot.
        figsize: Size of the figure.
        marker: Marker style for data points.
        rotation: Rotation angle for x-axis labels.
    """
    x_labels = [str(item[x_key]) for item in series]
    y_values = [item[y_key] for item in series]
    x_indices = list(range(len(series)))
    figure, ax = plt.subplots(figsize=figsize)
    ax.plot(x_indices, y_values, marker=marker)
    ax.set_xticks(x_indices)
    ax.set_xticklabels(x_labels, rotation=rotation, ha="right")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    plt.tight_layout()
    render_figure(figure)
