from typing import Dict, List, Union

import plotly.graph_objs as go
from mat3ra.made.material import Material
from mat3ra.made.tools.build.interface import StrainModes


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
