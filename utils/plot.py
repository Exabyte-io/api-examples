from typing import Dict, List, Union

import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objs as go
from ase.atoms import Atoms as ASEAtoms
from ase.optimize import BFGS, FIRE
from IPython.display import display
from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.rdf import RadialDistributionFunction
from mat3ra.made.tools.build.interface.enums import StrainModes
from plotly.subplots import make_subplots


def plot_strain_vs_atoms(interfaces: List[Material], settings: Dict[str, Union[str, int]]):
    """
    Plot the strain vs the number of atoms for each interface.
    Allows for visual inspection to find the most optimal interface for current task.
    Args:
        interfaces (List[Material]): The interfaces to plot.
        settings (Dict[str, Union[str, int]]): The settings for the plot.
    """

    data = []
    for index, interface in enumerate(interfaces):
        strain_percentage = interface.metadata[StrainModes.mean_abs_strain] * 100
        num_sites = len(interface.basis.coordinates.values)

        hover_text = f"Index: {index}<br>Strain: {strain_percentage:.2f}%<br>Atoms: {num_sites}"

        trace = go.Scatter(
            x=[strain_percentage],
            y=[num_sites],
            text=[hover_text],
            mode="markers",
            hoverinfo="text",
            name=f"Index: {index}",
        )
        data.append(trace)

    layout = go.Layout(
        xaxis=dict(title="Strain (%)", type=settings["X_SCALE"]),
        yaxis=dict(title="Number of atoms", type=settings["Y_SCALE"]),
        hovermode="closest",
        height=settings["HEIGHT"],
        legend_title_text="Interfaces Indices",
    )

    fig = go.Figure(data=data, layout=layout)
    fig.show()


def plot_twisted_interface_solutions(interfaces: List[Material]):
    """
    Plot the twisted interface solutions for number of atoms vs twist angle.
    Args:
        interfaces (List[Material]): The interfaces to plot.
    """
    data = []
    for i, interface in enumerate(interfaces):
        angle = interface.metadata.get("actual_twist_angle", 0)
        size = len(interface.basis.elements.ids)

        hover_text = f"Interface {i+1}<br>" f"Angle: {angle:.2f}°<br>" f"Atoms: {size}<br>"

        trace = go.Scatter(
            x=[angle], y=[size], text=[hover_text], mode="markers", hoverinfo="text", name=f"Interface {i+1}"
        )
        data.append(trace)

    layout = go.Layout(
        title="Twisted Interface Solutions",
        xaxis=dict(title="Twist Angle (°)"),
        yaxis=dict(title="Number of Atoms"),
        hovermode="closest",
    )

    fig = go.Figure(data=data, layout=layout)
    fig.show()


def create_realtime_plot():
    """
    Create a real-time plot for optimization progress.
    Returns:
        go.FigureWidget: The real-time plot.
    """
    fig = make_subplots(rows=1, cols=1, specs=[[{"type": "scatter"}]])
    scatter = go.Scatter(x=[], y=[], mode="lines+markers", name="Energy")
    fig.add_trace(scatter)
    fig.update_layout(title_text="Real-time Optimization Progress", xaxis_title="Step", yaxis_title="Energy (eV)")
    f = go.FigureWidget(fig)
    display(f)
    return f


def plot_update_callback(
    dynamic_object: Union[BFGS, FIRE],
    ase_interface: ASEAtoms,
    plotly_figure: go.FigureWidget,
    steps: List[int],
    energies: List[int],
):
    """
    Callback function for updating energies for steps in real-time.
    Args:
        dynamic_object: The ASE dynamics object.
        ase_interface: The ASE interface object.
        plotly_figure: The plotly figure widget.
        steps: The list of steps.
        energies: The list of energies.

    Returns:
        function: The update function that is attached to the dynamics object and called each step.
    """

    def update():
        step = dynamic_object.nsteps
        energy = ase_interface.get_total_energy()

        steps.append(step)
        energies.append(energy)

        print(f"Step: {step}, Energy: {energy:.4f} eV")
        with plotly_figure.batch_update():
            plotly_figure.data[0].x = steps
            plotly_figure.data[0].y = energies

    return update


def plot_rdf(material: Material, cutoff: float = 10.0, bin_size: float = 0.1):
    """
    Compute and plot the Radial Distribution Function (RDF) for a given material.

    Parameters:
    - material: The input material.
    - cutoff (float): Maximum distance for RDF calculation.
    - bin_size (float): Size of each bin in the histogram.

    Returns:
    - None
    """
    rdf = RadialDistributionFunction.from_material(material, cutoff=cutoff, bin_size=bin_size)

    # Plot the RDF
    plt.figure(figsize=(8, 5))
    plt.plot(rdf.bin_centers, rdf.rdf, label="Radial Distribution Function")
    plt.xlabel("Distance (Å)")
    plt.ylabel("g(r)")
    plt.title("Radial Distribution Function (RDF)")
    plt.legend()
    plt.grid()
    plt.show()


def plot_energy_landscape(xy_matrix, energy_matrix, optimal_position=None):
    """
    Create a 3D surface plot of the energy landscape.

    Args:
        xy_matrix (List[np.ndarray]): X and Y coordinate matrices
        energy_matrix (np.ndarray): Matrix of energy values
        optimal_position (tuple, optional): The optimal (x,y) position to highlight
    """
    x_vals, y_vals = xy_matrix
    X, Y = np.meshgrid(x_vals, y_vals)

    # Create the 3D surface plot
    fig = go.Figure(data=[go.Surface(x=X, y=Y, z=energy_matrix, colorscale="Viridis")])

    # Add optimal position marker if provided
    if optimal_position is not None:
        x_opt, y_opt = optimal_position[0], optimal_position[1]
        z_opt = np.min(energy_matrix)
        fig.add_trace(
            go.Scatter3d(
                x=[x_opt],
                y=[y_opt],
                z=[z_opt],
                mode="markers",
                marker=dict(size=8, color="red"),
                name="Optimal Position",
            )
        )

    fig.update_layout(
        title="Interface Energy Landscape",
        scene=dict(xaxis_title="X Position", yaxis_title="Y Position", zaxis_title="Energy"),
        width=800,
        height=800,
    )

    fig.show()


def plot_energy_heatmap(xy_matrix, energy_matrix, optimal_position=None):
    """
    Create a 2D heatmap of the energy landscape.

    Args:
        xy_matrix (List[np.ndarray]): X and Y coordinate matrices
        energy_matrix (np.ndarray): Matrix of energy values
        optimal_position (tuple, optional): The optimal (x,y) position to highlight
    """
    x_vals, y_vals = xy_matrix

    fig = go.Figure(
        data=go.Heatmap(x=x_vals, y=y_vals, z=energy_matrix, colorscale="Viridis", colorbar=dict(title="Energy"))
    )

    if optimal_position is not None:
        x_opt, y_opt = optimal_position[0], optimal_position[1]
        fig.add_trace(
            go.Scatter(
                x=[x_opt],
                y=[y_opt],
                mode="markers",
                marker=dict(size=12, color="red", symbol="x"),
                name="Optimal Position",
            )
        )

    fig.update_layout(
        title="Interface Energy Heatmap", xaxis_title="X Position", yaxis_title="Y Position", width=800, height=600
    )

    fig.show()
