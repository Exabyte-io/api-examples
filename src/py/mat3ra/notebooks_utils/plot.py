from matplotlib import pyplot as plt

from .ipython.plot._matplotlib import display_matplotlib_figure, plot_distribution_function, plot_series
from .ipython.plot._plotly import (
    create_realtime_plot,
    create_scatter_plot_2d,
    create_update_callback,
    plot_2d_heatmap,
    plot_3d_surface,
)
from .primitive.environment import is_pyodide_environment


def configure_matplotlib_renderer() -> None:
    if is_pyodide_environment():
        plt.switch_backend("Agg")


configure_matplotlib_renderer()

__all__ = [
    "display_matplotlib_figure",
    "plot_distribution_function",
    "plot_series",
    "create_scatter_plot_2d",
    "create_realtime_plot",
    "create_update_callback",
    "plot_2d_heatmap",
    "plot_3d_surface",
]
