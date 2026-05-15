from typing import Dict

from src.py.mat3ra.notebooks_utils.python.core.prompt import select_coordination_threshold_python

from .environment import is_pyodide_environment
from .jupyterlite.ui import select_coordination_threshold_emscripten
from .python.other.dataframe import dataframe_to_html


async def select_coordination_threshold(distribution: Dict[int, int], default_threshold: int) -> int:
    """
    Select the coordination threshold from the given distribution.
    Args:
        distribution:  The distribution of coordination numbers.
        default_threshold:  The default threshold value.
    Returns:
        int: The selected coordination threshold.
    """
    if is_pyodide_environment():
        return await select_coordination_threshold_emscripten(distribution, default_threshold)
    else:
        return select_coordination_threshold_python(distribution, default_threshold)


__all__ = ["dataframe_to_html"]
