import io
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objs as go
from IPython.display import Image, clear_output, display
from matplotlib.figure import Figure as MatplotlibFigure
from plotly.subplots import make_subplots

from .environment import is_pyodide_environment


def configure_matplotlib_renderer() -> None:
    if is_pyodide_environment():
        plt.switch_backend("Agg")


configure_matplotlib_renderer()


def display_matplotlib_figure(figure: MatplotlibFigure) -> None:
    buffer = io.BytesIO()
    figure.savefig(buffer, format="png")
    buffer.seek(0)
    display(Image(buffer.read()))
    plt.close(figure)


def render_figure(figure: Union[MatplotlibFigure, go.Figure, go.FigureWidget]) -> None:
    if is_pyodide_environment() and isinstance(figure, MatplotlibFigure):
        display_matplotlib_figure(figure)
        return
    if is_pyodide_environment():
        display(figure)
        return
    figure.show()


def scatter_plot_2d(
    x_values: List[float],
    y_values: List[float],
    hover_texts: List[str],
    settings: Dict[str, Any],
    trace_names: Optional[List[str]] = None,
) -> go.Figure:
    """
    Create a generic 2D scatter plot.

    Args:
        x_values: List of x-coordinates
        y_values: List of y-coordinates
        hover_texts: List of hover texts for each point
        settings: Plot settings including scales, height, and titles
        trace_names: Optional list of names for each trace
    """
    data = []
    for i in range(len(x_values)):
        trace = go.Scatter(
            x=[x_values[i]],
            y=[y_values[i]],
            text=[hover_texts[i]],
            mode="markers",
            hoverinfo="text",
            name=trace_names[i] if trace_names else f"Point {i}",
        )
        data.append(trace)

    layout = go.Layout(
        xaxis=dict(title=settings.get("x_title", "X"), type=settings.get("x_scale", "linear")),
        yaxis=dict(title=settings.get("y_title", "Y"), type=settings.get("y_scale", "linear")),
        hovermode="closest",
        height=settings.get("height", 600),
        title=settings.get("title", ""),
        legend_title_text=settings.get("legend_title", ""),
    )

    return go.Figure(data=data, layout=layout)


def create_realtime_plot(title: str = "Real-time Progress", x_label: str = "Step", y_label: str = "Value") -> go.Figure:
    """
    Create a real-time updating plot.
    """
    fig = make_subplots(rows=1, cols=1, specs=[[{"type": "scatter"}]])
    scatter = go.Scatter(x=[], y=[], mode="lines+markers", name="Progress")
    fig.add_trace(scatter)
    fig.update_layout(title_text=title, xaxis_title=x_label, yaxis_title=y_label)
    return fig


def create_update_callback(
    dynamic_object: Any,
    value_getter: Union[Callable, Any],
    figure: go.FigureWidget,
    steps: List[int],
    values: List[float],
    value_attr: Optional[str] = None,
    step_attr: str = "nsteps",
    print_format: str = "Step: {}, Value: {:.4f}",
):
    """
    Create a general update callback for real-time plotting.

    Args:
        dynamic_object: Object containing step information
        value_getter: Function to retrieve the measured value
        figure: Plotly figure to update
        steps: List to store step values
        values: List to store measured values
        step_attr: Attribute name for step count in dynamic_object
        value_attr: Optional attribute name for retrieving value from dynamic_object
        print_format: Format string for progress printing
    """

    def update():
        step = getattr(dynamic_object, step_attr, len(steps))

        if callable(value_getter):
            value = value_getter()
        elif value_attr:
            value = getattr(dynamic_object, value_attr, None)
        else:
            raise ValueError("Either value_getter (function) or value_attr (object attribute) must be provided.")

        steps.append(step)
        values.append(value)

        print(print_format.format(step, value))

        figure.data[0].x = steps
        figure.data[0].y = values

        # Update the plot by clearing and redrawing
        clear_output(wait=True)
        render_figure(figure)

    return update


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
    render_figure(figure)


def plot_3d_surface(
    x_matrix: np.ndarray,
    y_matrix: np.ndarray,
    z_matrix: np.ndarray,
    optimal_point: Optional[Tuple[float, float]] = None,
    title: str = "Surface Plot",
    labels: Optional[Dict[str, str]] = None,
) -> None:
    """
    Create a 3D surface plot with optional optimal point.

    Args:
        x_matrix: The x-axis matrix.
        y_matrix: The y-axis matrix.
        z_matrix: The z-axis matrix.
        optimal_point: The optimal point to highlight.
        title: The title of the plot.
        labels: The labels for the axes.
    """
    if labels is None:
        labels = {"x": "X", "y": "Y", "z": "Z"}

    fig = go.Figure(data=[go.Surface(x=x_matrix, y=y_matrix, z=z_matrix, colorscale="Viridis")])

    if optimal_point is not None:
        x_opt, y_opt = optimal_point
        z_opt = np.min(z_matrix)
        fig.add_trace(
            go.Scatter3d(
                x=[x_opt], y=[y_opt], z=[z_opt], mode="markers", marker=dict(size=8, color="red"), name="Optimal Point"
            )
        )

    fig.update_layout(
        title=title,
        scene=dict(xaxis_title=labels["x"], yaxis_title=labels["y"], zaxis_title=labels["z"]),
        width=800,
        height=800,
    )
    render_figure(fig)


def plot_2d_heatmap(
    x_values: np.ndarray,
    y_values: np.ndarray,
    z_matrix: np.ndarray,
    optimal_point: Optional[Tuple[float, float]] = None,
    title: str = "Heatmap",
    labels: Optional[Dict[str, str]] = None,
) -> None:
    """
    Create a 2D heatmap with optional optimal point.

    Args:
        x_values: The x-axis values.
        y_values: The y-axis values.
        z_matrix: The z-axis matrix.
        optimal_point: The optimal point to highlight.
        title: The title of the plot.
        labels: The labels for the axes.
    """
    if labels is None:
        labels = {"x": "X", "y": "Y", "z": "Z"}

    fig = go.Figure(
        data=go.Heatmap(x=x_values, y=y_values, z=z_matrix, colorscale="Viridis", colorbar=dict(title=labels["z"]))
    )

    if optimal_point is not None:
        x_opt, y_opt = optimal_point
        fig.add_trace(
            go.Scatter(
                x=[x_opt],
                y=[y_opt],
                mode="markers",
                marker=dict(size=12, color="red", symbol="x"),
                name="Optimal Point",
            )
        )

    fig.update_layout(title=title, xaxis_title=labels["x"], yaxis_title=labels["y"], width=800, height=600)
    render_figure(fig)
