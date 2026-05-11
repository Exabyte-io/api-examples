import json
import os
import sys
import uuid
from typing import Any, Dict, List, Optional, Union

from IPython.display import HTML, display
from pandas import DataFrame
from pandas.io.formats.style import Styler

from utils import settings


def create_prompt_text(array: List[Any], element_name: str = "element", prompt_head: Optional[str] = None) -> str:
    """
    Create a prompt text for selecting an element from an array.
    Args:
        array (List[Any]): The array (list) of elements to select from.
        element_name (str): The name of the element to be used in the prompt text (e.g., "transformation", "interface")
        prompt_head: The prompt text to be displayed before the list of elements.

    Returns:
        str: The prompt text.
    """
    prompt_head = prompt_head or f"Select {element_name} by index:\n"
    prompt_body = "\n".join(f"{i}: {t}" for i, t in enumerate(array))
    return prompt_head + prompt_body


def get_integer_from_input(selected_index_str: str, array: List[Any]) -> Union[int, None]:
    """
    Get an integer from the input string and check if it is a valid index for the given array.
    Args:
        selected_index_str (str): The input string to convert to an integer.
        array (List[Any]): The array (list) of elements to select from.

    Returns:
        int: The selected index if it is valid, otherwise None.
    """
    try:
        selected_index = int(selected_index_str)
    except ValueError:
        print("Invalid input. Please enter a valid integer.")
        return None

    if selected_index < 0 or selected_index >= len(array):
        print("Invalid index.")
        return None
    return selected_index


def ui_prompt_select_array_element_by_index(
    array: List[Any], element_name: str = "element", prompt_head: Optional[str] = None
) -> Any:
    """
    Prompt the user to select an element from an array by index and return the chosen element.
    Args:
        array (List[Any]): The array (list) of elements to select from.
        element_name (str): The name of the element to be used in the prompt text (e.g., "transformation", "interface")
        prompt_head: The prompt text to be displayed before the list of elements.

    Returns:
        Any: The selected element from the array.
    """
    prompt_text = create_prompt_text(array, element_name, prompt_head)
    selected_index_str = input(prompt_text)
    index = get_integer_from_input(selected_index_str, array)
    if index is None:
        return None
    result = array[index]
    print(f"Selected {element_name}: ", array[index])
    return result


async def ui_prompt_select_array_element_by_index_pyodide(
    array: List[Any], element_name: str = "element", prompt_head: Optional[str] = None
) -> Any:
    """
    Prompt the user to select an element from an array by index and return the chosen element.
    Used in Pyodide environment, needs to be awaited since the return of input() is type of PyodideFuture.
    Args:
        array (List[Any]): The array (list) of elements to select from.
        element_name (str): The name of the element to be used in the prompt text (e.g., "transformation", "interface")
        prompt_head: The prompt text to be displayed before the list of elements.

    Returns:
        Any: The selected element from the array.
    """
    prompt_text = create_prompt_text(array, element_name, prompt_head)
    # `input()` in Pyodide returns PyodideFuture, which is not compatible with `int` in regular Python environment
    selected_index_str = await input(prompt_text)  # type: ignore
    index = get_integer_from_input(selected_index_str, array)
    if index is None:
        return None
    result = array[index]
    print(f"Selected {element_name}: ", array[index])
    return result


async def select_coordination_threshold_emscripten(distribution: Dict[int, int], default_threshold: int) -> int:
    """
    Select the coordination threshold from the given distribution. Works in Pyodide environment.
    Args:
        distribution: The distribution of coordination numbers.
        default_threshold: The default threshold value.
    Returns:
        int: The selected coordination threshold.
    """
    coordination_threshold = default_threshold
    coordination_numbers = list(distribution.keys())
    prompt_text = f"\nCoordination numbers distribution: {distribution}" f"\nEnter coordination threshold value: "
    while True:
        try:
            value_str = await input(prompt_text)  # type: ignore
            value = int(value_str)
            if value in coordination_numbers:
                coordination_threshold = value
                break
            else:
                print(f"Invalid value. Please enter one of these coordination numbers: {coordination_numbers}")
                break
        except ValueError:
            print(f"Please enter a valid integer value from: {coordination_numbers}")
    return coordination_threshold


def select_coordination_threshold_python(distribution: Dict[int, int], default_threshold: int) -> int:
    """
    Select the coordination threshold from the given distribution. Works in regular Python environment.
    Args:
        distribution:  The distribution of coordination numbers.
        default_threshold:  The default threshold value.
    Returns:
        int: The selected coordination threshold.
    """
    coordination_threshold = default_threshold
    coordination_numbers = list(distribution.keys())
    prompt_text = f"\nCoordination numbers distribution: {distribution}" f"\nEnter coordination threshold value: "
    while True:
        try:
            value_str = input(prompt_text)
            value = int(value_str)
            if value in coordination_numbers:
                coordination_threshold = value
                break
            else:
                print(f"Invalid value. Please enter one of these coordination numbers: {coordination_numbers}")
                break
        except ValueError:
            print(f"Please enter a valid integer value from: {coordination_numbers}")
    return coordination_threshold


async def select_coordination_threshold(distribution: Dict[int, int], default_threshold: int) -> int:
    """
    Select the coordination threshold from the given distribution.
    Args:
        distribution:  The distribution of coordination numbers.
        default_threshold:  The default threshold value.
    Returns:
        int: The selected coordination threshold.
    """
    if sys.platform == "emscripten":
        return await select_coordination_threshold_emscripten(distribution, default_threshold)
    else:
        return select_coordination_threshold_python(distribution, default_threshold)


def dataframe_to_html(df: DataFrame, text_align: str = "center") -> Styler:
    """
    Converts Pandas dataframe to HTML.
    See https://pandas.pydata.org/pandas-docs/stable/style.html for more information about styling.

    Args:
        df (pd.DataFrame): Pandas dataframe.
        text_align (str): text align. Defaults to center.
    """
    styles = [
        dict(selector="th", props=[("text-align", text_align)]),
        dict(selector="td", props=[("text-align", text_align)]),
    ]
    return df.style.set_table_styles(styles)


def display_JSON(
    obj: Union[dict, list], interactive_viewer: bool = settings.use_interactive_JSON_viewer, level: int = 2
) -> None:
    """
    Displays JSON, either interactively or via a text dump to Stdout.

    The interactive viewer is based on https://github.com/mljar/mercury/blob/main/mercury/widgets/json.py.

    Args:
        obj (dict): Object to display as nicely-formatted JSON
        interactive_viewer (bool): Whether to use the interactive viewer or not
        level (int): The level to which the JSON should be expanded by default
    """
    if interactive_viewer:
        if isinstance(obj, (dict, list)):
            json_str = json.dumps(obj)
        else:
            json_str = obj

        id = str(uuid.uuid4())

        web_dir = os.path.join(os.path.dirname(__file__), "web")

        with open(os.path.join(web_dir, "renderjson.css")) as fp:
            css = fp.read()

        with open(os.path.join(web_dir, "renderjson.js")) as fp:
            js = fp.read()

        display(HTML(f'<style>{css}</style><div id="{id}"></div>'))
        display(
            HTML(
                f"<script>{js} "
                f"renderjson.set_show_to_level({str(level)}); "
                f'renderjson.set_icons("▸","▾"); '
                f'document.getElementById("{id}").appendChild(renderjson({json_str}))</script>'
            )
        )
    else:
        print(json.dumps(obj, indent=4))
