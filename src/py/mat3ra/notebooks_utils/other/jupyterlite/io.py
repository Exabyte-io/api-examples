import io
import os


async def write_to_file(file_name: str, file_content, mode: str = "wb"):
    """
    Write content to a file, handling both Python and Pyodide environments.
    In Python environment, the file is written to disk.
    In Pyodide environment, the file is written to `/drive/<folder-name>/<file-name>`.

    Args:
        file_name (str): The name of the file to write.
        file_content (str | bytes | io.StringIO | io.BytesIO): The content to write.
        mode (str): The mode to open the file in. Defaults to "wb" (write bytes).

    Returns:
        str: The absolute path of the saved file.
    """
    if isinstance(file_content, io.StringIO):
        file_content = file_content.getvalue().encode("utf-8")

    elif isinstance(file_content, io.BytesIO):
        file_content = file_content.getvalue()

    if "b" in mode and isinstance(file_content, str):
        file_content = file_content.encode("utf-8")

    with open(file_name, mode) as file:
        file.write(file_content)

    return os.path.abspath(file_name)
