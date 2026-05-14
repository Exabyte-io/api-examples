from typing import Any, Dict, List, Optional

from ..python.core.ui import create_prompt_text, get_integer_from_input


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
