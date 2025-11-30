import inspect
import json
import os
from typing import Any, Dict, List, Optional

from IPython.display import Javascript, display
from mat3ra.utils.jupyterlite.environment import ENVIRONMENT, EnvironmentsEnum
from mat3ra.utils.jupyterlite.logger import SeverityLevelEnum, log
from mat3ra.utils.jupyterlite.settings import UPLOADS_FOLDER


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


def set_data_python(key: str, value: Any):
    """
    Write data to the `uploads` folder in a JupyterLab environment.

    Args:
        key (str): The name under which data will be written.
        value (Any): The value to write to the `uploads` folder.
    """
    if not os.path.exists(UPLOADS_FOLDER):
        os.makedirs(UPLOADS_FOLDER)
    for item in value:
        safe_name = item["name"].replace("%", "pct").replace("/", ":")
        file_path = os.path.join(UPLOADS_FOLDER, f"{safe_name}.json")
        with open(file_path, "w") as file:
            json.dump(item, file)
        log(f"Data for {key} written to {file_path}")


def set_data(key: str, value: Any):
    """
    Switch between the two functions `set_data_pyodide` and `set_data_python` based on the environment.

    Args:
        key (str): The name under which data will be written or sent.
        value (Any): The value to write or send.
    """
    if ENVIRONMENT == EnvironmentsEnum.PYODIDE:
        set_data_pyodide(key, value)
    elif ENVIRONMENT == EnvironmentsEnum.PYTHON:
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


def get_data_python(key: str, globals_dict: Optional[Dict] = None):
    """
    Read data from the `uploads` folder in a JupyterLab environment.

    Args:
        key (str): The name under which data is expected to be received.
        globals_dict (dict, optional): A dictionary to store the received data. Defaults to None.
    """
    try:
        data_from_host = []
        index = 0
        for filename in sorted(os.listdir(UPLOADS_FOLDER)):
            if filename.endswith(".json"):
                with open(os.path.join(UPLOADS_FOLDER, filename), "r") as file:
                    data = json.load(file)
                name = os.path.splitext(filename)[0]
                log(f"{index}: {name}")
                index += 1
                data_from_host.append(data)
        if globals_dict is not None:
            globals_dict[key] = data_from_host
        return data_from_host
    except FileNotFoundError:
        print("No data found in the 'uploads' folder.")


def get_data(key: str, globals_dict: Optional[Dict] = None):
    """
    Switch between the two functions `get_data_pyodide` and `get_data_python` based on the environment.

    Args:
        key (str): The name under which data is expected to be received.
        globals_dict (dict, optional): A dictionary to store the received data. Defaults to None.
    """
    if ENVIRONMENT == EnvironmentsEnum.PYODIDE:
        get_data_pyodide(key, globals_dict)
    elif ENVIRONMENT == EnvironmentsEnum.PYTHON:
        get_data_python(key, globals_dict)


def get_materials(globals_dict: Optional[Dict] = None) -> List[Any]:
    """
    Retrieve materials from the environment and assign them to globals_dict["materials_in"].

    Args:
        globals_dict (dict, optional): The globals dictionary to populate.

    Returns:
        List[Material]: A list of Material objects.
    """
    from mat3ra.made.material import Material
    from mat3ra.made.tools.build_components import MaterialWithBuildMetadata

    if globals_dict is None:
        # Get the globals of the caller for correct variable assignment during the execution of data_bridge extension
        frame = inspect.currentframe()
        try:
            caller_frame = frame.f_back  # type: ignore
            caller_globals = caller_frame.f_globals  # type: ignore
            globals_dict = caller_globals
        finally:
            del frame  # Avoid reference cycles
    get_data("materials_in", globals_dict)

    if "materials_in" in globals_dict and globals_dict["materials_in"]:
        materials = []
        for item in globals_dict["materials_in"]:
            try:
                materials.append(MaterialWithBuildMetadata.create(item))
            except Exception:
                materials.append(Material.create(item))
        log(f"Retrieved {len(materials)} materials.")
        return materials
    else:
        # Fallback to load materials from the UPLOADS_FOLDER if launched outside of Materials Designer
        log(f"No input materials found. Loading from the {UPLOADS_FOLDER} folder.")
        return load_materials_from_folder()


def set_materials(materials: List[Any]):
    """
    Serialize and send a list of Material objects to the environment.

    Args:
        materials (List[Material]): The list of Material objects to send.
    """
    from mat3ra.utils.array import convert_to_array_if_not

    materials = convert_to_array_if_not(materials)
    materials_data = [json.loads(material.to_json()) for material in materials]
    set_data("materials", materials_data)


def load_materials_from_folder(folder_path: Optional[str] = None, verbose: bool = True) -> List[Any]:
    """
    Load materials from the specified folder or from the UPLOADS_FOLDER by default.

    Args:
        folder_path (Optional[str]): The path to the folder containing material files.
                                     If not provided, defaults to the UPLOADS_FOLDER.
        verbose (bool): Whether to log verbose messages.

    Returns:
        List[Material]: A list of Material objects loaded from the folder.
    """
    from mat3ra.made.material import Material
    from mat3ra.made.tools.build_components import MaterialWithBuildMetadata

    folder_path = folder_path or UPLOADS_FOLDER

    if not os.path.exists(folder_path):
        log(f"Folder '{folder_path}' does not exist.", SeverityLevelEnum.ERROR, force_verbose=verbose)
        return []

    data_from_host = []
    try:
        index = 0
        for filename in sorted(os.listdir(folder_path)):
            if filename.endswith(".json"):
                with open(os.path.join(folder_path, filename), "r") as file:
                    data = json.load(file)
                name = os.path.splitext(filename)[0]
                log(f"{index}: {name}", SeverityLevelEnum.INFO, force_verbose=verbose)
                index += 1
                data_from_host.append(data)
    except FileNotFoundError:
        log(f"No data found in the '{folder_path}' folder.", SeverityLevelEnum.ERROR, force_verbose=verbose)
        return []
    try:
        materials = [MaterialWithBuildMetadata.create(item) for item in data_from_host]
    except Exception:
        materials = [Material.create(item) for item in data_from_host]

    if materials:
        log(
            f"Successfully loaded {len(materials)} materials from folder '{folder_path}'",
            SeverityLevelEnum.INFO,
            force_verbose=verbose,
        )
    else:
        log(f"No materials found in folder '{folder_path}'", SeverityLevelEnum.WARNING, force_verbose=verbose)

    return materials


def load_material_from_folder(folder_path: str, name: str, verbose: bool = True) -> Optional[Any]:
    """
    Load a single material from the specified folder by matching a substring of the name.

    Args:
        folder_path (str): The path to the folder containing material files.
        name (str): The substring of the name of the material to load.
        verbose (bool): Whether to log verbose messages.

    Returns:
        Optional[Material]: The first Material object that contains the name substring, or None if not found.
    """
    # Reuse the existing function to load all materials from the folder
    materials = load_materials_from_folder(folder_path, verbose=verbose)
    for material in materials:
        if name.lower() in material.name.lower():
            log(
                f"Found: '{material.name}'",
                SeverityLevelEnum.INFO,
                force_verbose=verbose,
            )
            return material

    log(
        f"No material containing '{name}' found in folder '{folder_path}'.",
        SeverityLevelEnum.WARNING,
        force_verbose=verbose,
    )
    return None


def write_materials_to_folder(materials: List[Any], folder_path: Optional[str] = None, verbose: bool = True):
    """
    Write materials to the specified folder or to the UPLOADS_FOLDER by default.

    Args:
        materials (List[Material]): The list of Material objects to write to the folder.
        folder_path (Optional[str]): The path to the folder where the materials will be written.
                                     If not provided, defaults to the UPLOADS_FOLDER.
        verbose (bool): Whether to log verbose messages.
    """
    from mat3ra.utils.array import convert_to_array_if_not

    folder_path = folder_path or UPLOADS_FOLDER
    materials = convert_to_array_if_not(materials)

    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
        if verbose:
            log(f"Created folder '{folder_path}'.", SeverityLevelEnum.INFO)

    for material in materials:
        safe_name = material.name.replace("%", "pct").replace("/", ":")
        file_path = os.path.join(folder_path, f"{safe_name}.json")
        with open(file_path, "w") as file:
            json.dump(material.to_dict(), file)
        log(f"Material '{material.name}' written to '{file_path}'", SeverityLevelEnum.INFO, force_verbose=verbose)


def download_content_to_file(content: Any, filename: str):
    """
    Download content to a file with the given filename.

    Args:
        content (Any): The content to download.
        filename (str): The name of the file to download.
    """
    from mat3ra.made.material import Material
    from mat3ra.made.tools.build_components import MaterialWithBuildMetadata

    if isinstance(content, dict):
        content = json.dumps(content, indent=4)

    if isinstance(content, (Material, MaterialWithBuildMetadata)):
        content = content.to_json()
        content = json.dumps(content, indent=4)

    js_code = f"""
    var content = `{content}`;
    var filename = `{filename}`;
    var blob = new Blob([content], {{ type: 'application/json' }});
    var link = document.createElement('a');
    link.href = window.URL.createObjectURL(blob);
    link.download = filename;
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    """
    display(Javascript(js_code))
