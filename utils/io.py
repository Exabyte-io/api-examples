import sys
from typing import Any, List, Optional, Union


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


async def select_coordination_threshold_emscripten(coordination_numbers, default_threshold):
    coordination_threshold = default_threshold
    prompt_text = f"\nCoordination numbers: {coordination_numbers}" f"\nEnter coordination threshold value: "
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


def select_coordination_threshold_python(coordination_numbers, default_threshold):
    coordination_threshold = default_threshold
    prompt_text = f"\nCoordination numbers: {coordination_numbers}" f"\nEnter coordination threshold value: "
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


async def select_coordination_threshold(coordination_numbers, default_threshold):
    if sys.platform == "emscripten":
        return await select_coordination_threshold_emscripten(coordination_numbers, default_threshold)
    else:
        return select_coordination_threshold_python(coordination_numbers, default_threshold)
