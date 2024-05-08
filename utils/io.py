from typing import Any, Awaitable, Callable, List, Optional


def generate_prompt_text(array: List[Any]) -> str:
    return "\n".join(f"{i}: {t}" for i, t in enumerate(array))


async def async_input(prompt: str) -> str:
    return input(prompt)


def select_index(
    input_func: Callable[[str], Awaitable[str] or str],
    array: List[Any],
    element_name: str,
    prompt_head: Optional[str] = None,
) -> Any:
    prompt_head = prompt_head or f"Select {element_name} by index:\n"
    prompt_body = generate_prompt_text(array)
    selected_index_str = input_func(prompt_head + prompt_body)
    try:
        selected_index = int(selected_index_str)
    except ValueError:
        print("Invalid input. Please enter a valid integer.")
        return select_index(input_func, array, element_name, prompt_head)

    if selected_index < 0 or selected_index >= len(array):
        print("Invalid index.")
        return select_index(input_func, array, element_name, prompt_head)

    print(f"Selected {element_name}: ", array[selected_index])
    return array[selected_index]


def user_select_array_element_by_index(
    array: List[Any], element_name: str = "element", prompt_head: Optional[str] = None
) -> Any:
    return select_index(input, array, element_name, prompt_head)


async def user_select_array_element_by_index_async(
    array: List[Any], element_name: str = "element", prompt_head: Optional[str] = None
) -> Any:
    return await select_index(async_input, array, element_name, prompt_head)
