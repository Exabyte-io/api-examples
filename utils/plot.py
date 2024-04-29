from collections import defaultdict

import plotly.graph_objs as go


def plot_strain_vs_atoms(sorted_interfaces, settings):
    terminations = sorted_interfaces.terminations
    # Create a mapping from termination to its index
    termination_to_index = {termination: i for i, termination in enumerate(terminations)}

    grouped_interfaces = defaultdict(list)
    for termination, interfaces in sorted_interfaces.items():
        for index, interface_data in enumerate(interfaces):
            strain_percentage = interface_data["mean_abs_strain"] * 100
            num_sites = interface_data["interface"].num_sites
            key = (strain_percentage, num_sites)
            grouped_interfaces[key].append((index, termination))

    data = []
    for (strain, num_sites), indices_and_terminations in grouped_interfaces.items():
        termination_indices = defaultdict(list)
        for index, termination in indices_and_terminations:
            termination_indices[termination].append(index)
        all_indices = [index for indices in termination_indices.values() for index in indices]
        index_range = f"{min(all_indices)}-{max(all_indices)}" if len(all_indices) > 1 else str(min(all_indices))

        hover_text = "<br>-----<br>".join(
            f"Termination: {termination}<br>Termination index: {termination_to_index[termination]}"
            f"<br>Interfaces Index Range: {index_range}<br>Strain: {strain:.2f}%<br>Atoms: {num_sites}"
            for termination, indices in termination_indices.items()
        )
        trace = go.Scatter(
            x=[strain],
            y=[num_sites],
            text=[hover_text],
            mode="markers",
            hoverinfo="text",
            name=f"Indices: {index_range}",
        )
        data.append(trace)

    layout = go.Layout(
        xaxis=dict(title="Strain (%)", type=settings["X_SCALE"]),
        yaxis=dict(title="Number of atoms", type=settings["Y_SCALE"]),
        hovermode="closest",
        height=settings["HEIGHT"],
        legend_title_text="Interfaces Index Range",
    )
    fig = go.Figure(data=data, layout=layout)
    fig.show()
