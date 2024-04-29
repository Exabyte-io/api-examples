import plotly.graph_objs as go


def plot_strain_vs_atoms(interface_data_holder, settings):
    sorted_interfaces = interface_data_holder.interfaces

    data = []
    for termination, interfaces in sorted_interfaces.items():
        for index, interface in enumerate(interfaces):
            strain_percentage = interface.get_mean_abs_strain() * 100
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
        legend_title_text="Interfaces Index Range",
    )

    fig = go.Figure(data=data, layout=layout)
    fig.show()
