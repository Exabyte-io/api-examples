import io
import json
import os
from typing import Any, Dict, Optional

from IPython.display import Javascript, display

from ..logger import log
from ..python.io import set_data_python


def set_data_pyodide(key: str, value: Any):
    """
    Take a Python object, serialize it to JSON, and send it to the host environment
    through a JavaScript function defined in the JupyterLite extension `data_bridge`.

    Args:
        key (str): The name under which data will be sent.
        value (Any): The value to send to the host environment.
    """
    serialized_data = json.dumps({key: value})
    js_code = f"""
      (function() {{
          if (window.sendDataToHost) {{
              window.sendDataToHost({serialized_data});
              console.log('Data sent to host:', {serialized_data});
          }} else {{
              console.error('sendDataToHost function is not defined on the window object.');
          }}
      }})();
      """
    display(Javascript(js_code))
    log(f"Data for {key} sent to host.")
    set_data_python(key, value)


def get_data_pyodide(key: str, globals_dict: Optional[Dict] = None):
    """
    Load data from the host environment into globals()[key] variable.

    Args:
        key (str): Global variable name to store the received data.
        globals_dict (dict, optional): globals() dictionary of the current scope.
    """
    if globals_dict is not None:
        globals_dict[key] = globals_dict.get("data_from_host", None)


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
