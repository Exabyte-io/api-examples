from matplotlib import pyplot as plt

from .environment import is_pyodide_environment


def configure_matplotlib_renderer() -> None:
    if is_pyodide_environment():
        plt.switch_backend("Agg")
