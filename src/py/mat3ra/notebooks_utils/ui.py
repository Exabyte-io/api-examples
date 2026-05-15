from typing import Dict

from .core.prompt import select_coordination_threshold_emscripten
from .ipython.ui import dataframe_to_html, display_JSON
from .primitive.environment import is_pyodide_environment
from .primitive.prompt import select_coordination_threshold_python


async def select_coordination_threshold(distribution: Dict[int, int], default_threshold: int) -> int:
    """
    Select the coordination threshold from the given distribution.

    Args:
        distribution: The distribution of coordination numbers.
        default_threshold: The default threshold value.

    Returns:
        int: The selected coordination threshold.
    """
    if is_pyodide_environment():
        return await select_coordination_threshold_emscripten(distribution, default_threshold)
    else:
        return select_coordination_threshold_python(distribution, default_threshold)


__all__ = ["dataframe_to_html", "display_JSON"]
