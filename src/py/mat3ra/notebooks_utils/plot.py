from .ipython.plot._matplotlib import (
    configure_matplotlib_renderer,
    display_matplotlib_figure,
    plot_distribution_function,
    plot_series,
)
from .ipython.plot._plotly import (
    create_realtime_plot,
    create_scatter_plot_2d,
    create_update_callback,
    plot_2d_heatmap,
    plot_3d_surface,
)

configure_matplotlib_renderer()

__all__ = [
    "configure_matplotlib_renderer",
    "display_matplotlib_figure",
    "plot_distribution_function",
    "plot_series",
    "create_scatter_plot_2d",
    "create_realtime_plot",
    "create_update_callback",
    "plot_2d_heatmap",
    "plot_3d_surface",
]
