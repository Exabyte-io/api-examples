from typing import Dict, List, Union

import plotly.graph_objs as go
from IPython.display import display
from mat3ra.made.material import Material
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
        num_sites = len(interface.basis["coordinates"])

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


def create_realtime_plot():
    fig = make_subplots(rows=1, cols=1, specs=[[{"type": "scatter"}]])
    scatter = go.Scatter(x=[], y=[], mode="lines+markers", name="Energy")
    fig.add_trace(scatter)
    fig.update_layout(title_text="Real-time Optimization Progress", xaxis_title="Step", yaxis_title="Energy (eV)")
    f = go.FigureWidget(fig)
    display(f)
    return f


def update_plot(f, steps, energies):
    with f.batch_update():
        f.data[0].x = steps
        f.data[0].y = energies


def plot_update_callback(dyn, ase_interface, fig, steps, energies):
    def update():
        step = dyn.nsteps
        energy = ase_interface.get_total_energy()

        steps.append(step)
        energies.append(energy)

        print(f"Step: {step}, Energy: {energy:.4f} eV")
        update_plot(fig, steps, energies)

    return update
