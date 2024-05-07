from typing import Any, List, Optional


def user_select_object_by_index(
    objects: List[Any], object_name: str = "object", prompt_head: Optional[str] = None
) -> Any:
    prompt_head = prompt_head or f"Select {object_name} by index:\n"
    prompt_body = "\n".join(f"{i}: {t}" for i, t in enumerate(objects))

    selected_index_str = input(prompt_head + prompt_body)
    try:
        selected_index = int(selected_index_str)
    except ValueError:
        print("Invalid input. Please enter a valid integer.")
        return user_select_object_by_index(objects, object_name, prompt_head)

    assert isinstance(selected_index, int)

    if selected_index < 0 or selected_index >= len(objects):
        print("Invalid index.")
        return user_select_object_by_index(objects, object_name, prompt_head)

    print(f"Selected {object_name}: ", objects[selected_index])
    return objects[selected_index]
