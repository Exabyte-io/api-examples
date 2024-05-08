from typing import Any, List, Optional


def create_prompt_text(array: List[Any], element_name: str = "element", prompt_head: Optional[str] = None) -> str:
    prompt_head = prompt_head or f"Select {element_name} by index:\n"
    prompt_body = "\n".join(f"{i}: {t}" for i, t in enumerate(array))
    return prompt_head + prompt_body


def process_selected_index(selected_index_str: str, array: List[Any], element_name: str) -> Any:
    try:
        selected_index = int(selected_index_str)
    except ValueError:
        print("Invalid input. Please enter a valid integer.")
        return None  # Return None to indicate failure

    if selected_index < 0 or selected_index >= len(array):
        print("Invalid index.")
        return None  # Return None to indicate failure

    print(f"Selected {element_name}: ", array[selected_index])
    return array[selected_index]


def user_select_array_element_by_index(
    array: List[Any], element_name: str = "element", prompt_head: Optional[str] = None
) -> Any:
    prompt_text = create_prompt_text(array, element_name, prompt_head)
    selected_index_str = input(prompt_text)
    result = process_selected_index(selected_index_str, array, element_name)
    if result is None:
        return user_select_array_element_by_index(array, element_name, prompt_head)
    return result


async def user_select_array_element_by_index_pyodide(
    array: List[Any], element_name: str = "element", prompt_head: Optional[str] = None
) -> Any:
    prompt_text = create_prompt_text(array, element_name, prompt_head)
    # `input()` in Pyodide returns PyodideFuture, which is not compatible with `int` in regular Python environment
    selected_index_str = await input(prompt_text)  # type: ignore
    result = process_selected_index(selected_index_str, array, element_name)
    if result is None:
        return await user_select_array_element_by_index_pyodide(array, element_name, prompt_head)
    return result
