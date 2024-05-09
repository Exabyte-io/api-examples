from typing import Any, List, Optional


def create_prompt_text(array: List[Any], element_name: str = "element", prompt_head: Optional[str] = None) -> str:
    prompt_head = prompt_head or f"Select {element_name} by index:\n"
    prompt_body = "\n".join(f"{i}: {t}" for i, t in enumerate(array))
    return prompt_head + prompt_body


def get_integer_from_input(selected_index_str: str, array: List[Any]) -> Any:
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
    prompt_text = create_prompt_text(array, element_name, prompt_head)
    selected_index_str = input(prompt_text)
    index = get_integer_from_input(selected_index_str, array)
    if index is None:
        return ui_prompt_select_array_element_by_index(array, element_name, prompt_head)
    result = array[index]
    print(f"Selected {element_name}: ", array[index])
    return result


async def ui_prompt_select_array_element_by_index_pyodide(
    array: List[Any], element_name: str = "element", prompt_head: Optional[str] = None
) -> Any:
    prompt_text = create_prompt_text(array, element_name, prompt_head)
    # `input()` in Pyodide returns PyodideFuture, which is not compatible with `int` in regular Python environment
    selected_index_str = await input(prompt_text)  # type: ignore
    index = get_integer_from_input(selected_index_str, array)
    if index is None:
        return await ui_prompt_select_array_element_by_index_pyodide(array, element_name, prompt_head)
    result = array[index]
    print(f"Selected {element_name}: ", array[index])
    return result
