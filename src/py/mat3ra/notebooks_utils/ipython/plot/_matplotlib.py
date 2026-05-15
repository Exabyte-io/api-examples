import io
from typing import Dict, List, Tuple

import numpy as np
from IPython.display import Image, display
from matplotlib import pyplot as plt
from matplotlib.figure import Figure as MatplotlibFigure


def display_matplotlib_figure(figure: MatplotlibFigure) -> None:
    buffer = io.BytesIO()
    figure.savefig(buffer, format="png")
    buffer.seek(0)
    display(Image(buffer.read()))
    plt.close(figure)


def plot_distribution_function(
    bin_centers: np.ndarray,
    distribution: np.ndarray,
    xlabel: str = "Distance",
    ylabel: str = "g(r)",
    title: str = "Distribution Function",
    figsize: Tuple[int, int] = (8, 5),
) -> None:
    """
    Plot a generic distribution function.

    Args:
        bin_centers: The bin centers.
        distribution: The distribution values.
        xlabel: The x-axis label.
        ylabel: The y-axis label.
        title: The title of the plot.
        figsize: The size of the figure.
    """
    figure = plt.figure(figsize=figsize)
    plt.plot(bin_centers, distribution, label=title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.grid()
    display_matplotlib_figure(figure)


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
    display_matplotlib_figure(figure)
