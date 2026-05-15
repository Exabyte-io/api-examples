from typing import Dict


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
