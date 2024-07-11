import json
import os
from enum import Enum
from typing import Any, Dict, Optional

from IPython.display import Javascript, display

UPLOADS_FOLDER = "uploads"


class EnvironmentEnum(Enum):
    PYODIDE = "pyodide"
    PYTHON = "python"


# Environment detection
# default value for env.HOME from https://pyodide.org/en/stable/usage/api/js-api.html
ENVIRONMENT = EnvironmentEnum.PYODIDE if os.environ.get("HOME") == "/home/pyodide" else EnvironmentEnum.PYTHON

if ENVIRONMENT == EnvironmentEnum.PYODIDE:
    import micropip

if ENVIRONMENT == EnvironmentEnum.PYTHON:
    import subprocess
    import sys


async def install_package_pyodide(pkg: str, verbose=True):
    """
    Install a package in a Pyodide environment.
    Args:
        pkg (string): The name of the package to install.
        verbose (bool): Whether to print the name of the installed package.
    """
    is_url = pkg.startswith("http://") or pkg.startswith("https://") or pkg.startswith("emfs:/")
    are_dependencies_installed = not is_url
    await micropip.install(pkg, deps=are_dependencies_installed)
    pkg_name = pkg.split("/")[-1].split("-")[0] if is_url else pkg.split("==")[0]
    if verbose:
        print(f"Installed {pkg_name}")


def install_package_python(pkg: str, verbose=True):
    """
    Install a package in a standard Python environment.
    Args:
        pkg (string): The name of the package to install.
        verbose (bool): Whether to print the name of the installed package.
    """
    subprocess.check_call([sys.executable, "-m", "pip", "install", pkg])
    if verbose:
        print(f"Installed {pkg}")


async def install_packages(notebook_name: str, requirements_path="config.yml", verbose=True):
    """
    Install the packages listed in the requirements file for the notebook with the given name.
    Args:
        notebook_name (string): The name of the notebook for which to install packages.
        requirements_path (string): The path to the requirements file.
        verbose (bool): Whether to print the names of the installed packages and status of installation.
    """
    if ENVIRONMENT == EnvironmentEnum.PYODIDE:
        await micropip.install("pyyaml")
        # PyYAML has to be installed before being imported in Pyodide and can't appear at the top of the file
    import yaml

    base_path = os.getcwd()
    if requirements_path is None:
        requirements_file = os.path.normpath(os.path.join(base_path, "./config.yml"))
    else:
        requirements_file = os.path.normpath(os.path.join(base_path, requirements_path))

    with open(requirements_file, "r") as f:
        requirements = yaml.safe_load(f)

    # Hash the requirements to avoid re-installing packages
    requirements_hash = str(hash(json.dumps(requirements)))
    if os.environ.get("requirements_hash") != requirements_hash:
        packages_default_common = requirements.get("default", {}).get("packages_common", []) or []
        packages_default_environment_specific = (
            requirements.get("default", {}).get(f"packages_{ENVIRONMENT.value}", []) or []
        )

        notebook_requirements = next(
            (cfg for cfg in requirements.get("notebooks", []) if cfg.get("name") == notebook_name), None
        )
        if notebook_requirements:
            packages_notebook_common = notebook_requirements.get("packages_common", []) or []
            packages_notebook_environment_specific = (
                notebook_requirements.get(f"packages_{ENVIRONMENT.value}", []) or []
            )
        else:
            packages_notebook_common = []
            packages_notebook_environment_specific = []

        # Note: environment specific packages have to be installed first,
        # because in Pyodide common packages might depend on them
        packages = [
            *packages_default_environment_specific,
            *packages_notebook_environment_specific,
            *packages_default_common,
            *packages_notebook_common,
        ]

        for pkg in packages:
            if ENVIRONMENT == EnvironmentEnum.PYODIDE:
                await install_package_pyodide(pkg, verbose)
            elif ENVIRONMENT == EnvironmentEnum.PYTHON:
                install_package_python(pkg, verbose)

        if verbose:
            print("Packages installed successfully.")
        os.environ["requirements_hash"] = requirements_hash


def set_data_pyodide(key: str, value: Any):
    """
    Take a Python object, serialize it to JSON, and send it to the host environment
    through a JavaScript function defined in the JupyterLite extension `data_bridge`.

    Args:
        key (string): The name under which data will be sent.
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
    print(f"Status: {key} sent to host.")
    set_data_python(key, value)


def set_data_python(key: str, value: Any):
    """
    Write data to the `uploads` folder in a JupyterLab environment.
    Args:
        key (string): The name under which data will be written.
        value (Any): The value to write to the `uploads` folder.
    """
    if not os.path.exists(UPLOADS_FOLDER):
        os.makedirs(UPLOADS_FOLDER)
    for item in value:
        safe_name = item["name"].replace("%", "pct").replace("/", ":")
        file_path = os.path.join(UPLOADS_FOLDER, f"{safe_name}.json")
        with open(file_path, "w") as file:
            json.dump(item, file)
        print(f"Data for {key} written to {file_path}")


def set_data(key: str, value: Any):
    """
    Switch between the two functions `set_data_pyodide` and `set_data_python` based on the environment.
    Args:
        key (string): The name under which data will be written or sent.
        value (Any): The value to write or send.
    """
    if ENVIRONMENT == EnvironmentEnum.PYODIDE:
        set_data_pyodide(key, value)
    elif ENVIRONMENT == EnvironmentEnum.PYTHON:
        set_data_python(key, value)


def get_data_pyodide(key: str, globals_dict: Optional[Dict] = None):
    """
    Load data from the host environment into globals()[key] variable.
    Args:
        key (string): global variable name to store the received data.
        globals_dict (dict): globals() dictionary of the current scope.
    """
    if globals_dict is not None:
        globals_dict[key] = globals_dict["data_from_host"]


def get_data_python(key: str, globals_dict: Optional[Dict] = None):
    """
    Read data from the `uploads` folder in a JupyterLab environment.
    Args:
        key (string): The name under which data is expected to be received.
        globals_dict (dict): A dictionary to store the received data. Defaults to None.
    """
    try:
        data_from_host = []
        for filename in sorted(os.listdir(UPLOADS_FOLDER)):
            if filename.endswith(".json"):
                with open(os.path.join(UPLOADS_FOLDER, filename), "r") as file:
                    data = json.load(file)
                name = os.path.splitext(filename)[0]
                print(f"Data from {name} has been read successfully.")
                data_from_host.append(data)
        if globals_dict is not None:
            globals_dict[key] = data_from_host
    except FileNotFoundError:
        print("No data found in the 'uploads' folder.")


def get_data(key: str, globals_dict: Optional[Dict] = None):
    """
    Switch between the two functions `get_data_pyodide` and `get_data_python` based on the environment.
    Args:
        key (string): The name under which data is expected to be received.
        globals_dict (dict): A dictionary to store the received data. Defaults to None.
    """
    if ENVIRONMENT == EnvironmentEnum.PYODIDE:
        get_data_pyodide(key, globals_dict)
    elif ENVIRONMENT == EnvironmentEnum.PYTHON:
        get_data_python(key, globals_dict)
