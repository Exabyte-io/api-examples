from typing import Any, List, Optional


def user_select_array_element_by_index(
    array: List[Any], element_name: str = "element", prompt_head: Optional[str] = None
) -> Any:
    prompt_head = prompt_head or f"Select {element_name} by index:\n"
    prompt_body = "\n".join(f"{i}: {t}" for i, t in enumerate(array))

    selected_index_str = input(prompt_head + prompt_body)
    try:
        selected_index = int(selected_index_str)
    except ValueError:
        print("Invalid input. Please enter a valid integer.")
        return user_select_array_element_by_index(array, element_name, prompt_head)

    assert isinstance(selected_index, int)

    if selected_index < 0 or selected_index >= len(array):
        print("Invalid index.")
        return user_select_array_element_by_index(array, element_name, prompt_head)

    print(f"Selected {element_name}: ", array[selected_index])
    return array[selected_index]


async def user_select_array_element_by_index_async(
    array: List[Any], element_name: str = "element", prompt_head: Optional[str] = None
) -> Any:
    prompt_head = prompt_head or f"Select {element_name} by index:\n"
    prompt_body = "\n".join(f"{i}: {t}" for i, t in enumerate(array))

    # `input` in Pyodide returns PyodideFuture, which is not compatible with `int` in regular Python environment
    selected_index_str = await input(prompt_head + prompt_body)  # type: ignore
    try:
        selected_index = int(selected_index_str)
    except ValueError:
        print("Invalid input. Please enter a valid integer.")
        return await user_select_array_element_by_index(array, element_name, prompt_head)

    assert isinstance(selected_index, int)

    if selected_index < 0 or selected_index >= len(array):
        print("Invalid index.")
        return await user_select_array_element_by_index(array, element_name, prompt_head)

    print(f"Selected {element_name}: ", array[selected_index])
    return array[selected_index]
