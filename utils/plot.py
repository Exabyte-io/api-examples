from typing import Dict, Union

import plotly.graph_objs as go
from mat3ra.made.tools.build.interface import InterfaceDataHolder, StrainModes


def plot_strain_vs_atoms(interface_data_holder: InterfaceDataHolder, settings: Dict[str, Union[str, int]]):
    """
    Plot the strain vs the number of atoms for each interface.
    Allows for visual inspection to find the most optimal interface for current task.
    Args:
        interface_data_holder (InterfaceDataHolder): The interface data holder object.
        settings (Dict[str, Union[str, int]]): The settings for the plot.
    """
    sorted_interfaces = interface_data_holder.get_interfaces_for_termination(0)

    data = []
    for index, interface in enumerate(sorted_interfaces):
        strain_percentage = interface.interface_properties[StrainModes.mean_abs_strain] * 100
        num_sites = interface.num_sites

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
