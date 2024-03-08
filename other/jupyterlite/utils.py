from IPython.display import display, Javascript
import json
import yaml
import os
from enum import Enum

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


async def install_package_pyodide(pkg, verbose=True):
    """
    Install a package in a Pyodide environment.
    Args:
        pkg (string): The name of the package to install.
        verbose (bool): Whether to print the name of the installed package.
    """
    is_url = pkg.startswith("http://") or pkg.startswith("https://")
    are_dependencies_installed = not is_url
    await micropip.install(pkg, deps=are_dependencies_installed)
    pkg_name = pkg.split("/")[-1].split("-")[0] if is_url else pkg.split("==")[0]
    if verbose:
        print(f"Installed {pkg_name}")


def install_package_python(pkg, verbose=True):
    """
    Install a package in a standard Python environment.
    Args:
        pkg (string): The name of the package to install.
        verbose (bool): Whether to print the name of the installed package.
    """
    subprocess.check_call([sys.executable, "-m", "pip", "install", pkg])
    if verbose:
        print(f"Installed {pkg}")


# Set the install_package function based on the environment
install_package = install_package_pyodide if ENVIRONMENT == EnvironmentEnum.PYODIDE else install_package_python


async def install_packages(notebook_name, requirements_path="config.yml", verbose=True):
    """
    Install the packages listed in the requirements file for the notebook with the given name.
    Args:
        notebook_name (string): The name of the notebook for which to install packages.
        requirements_path (string): The path to the requirements file.
        verbose (bool): Whether to print the names of the installed packages and status of installation.
    """
    with open(requirements_path, "r") as f:
        requirements = yaml.safe_load(f)

    # Hash the requirements to avoid re-installing packages
    requirements_hash = str(hash(json.dumps(requirements)))
    if os.environ.get("requirements_hash") != requirements_hash:
        packages = requirements.get("default", {}).get("packages_common", [])
        packages += requirements.get("default", {}).get(f"packages_{ENVIRONMENT.value}", [])

        notebook_requirements = next(
            (cfg for cfg in requirements.get("notebooks", []) if cfg.get("name") == notebook_name), None
        )
        if notebook_requirements:
            packages += notebook_requirements.get("packages_common", [])
            packages += notebook_requirements.get(f"packages_{ENVIRONMENT.value}", [])
        else:
            raise ValueError(f"No packages found for notebook {notebook_name}")

        for pkg in packages:
            await install_package(pkg, verbose)

        if verbose:
            print("Packages installed successfully.")
        os.environ["requirements_hash"] = requirements_hash


def set_data_pyodide(key, value):
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


def set_data_python(key, value):
    """
    Write data to the `uploads` folder in a JupyterLab environment.
    Args:
        key (string): The name under which data will be written.
        value (Any): The value to write to the `uploads` folder.
    """
    if not os.path.exists(UPLOADS_FOLDER):
        os.makedirs(UPLOADS_FOLDER)
    file_path = os.path.join(UPLOADS_FOLDER, f"{key}.json")
    with open(file_path, "w") as file:
        json.dump(value, file)
    print(f"Data for {key} written to {file_path}")


def set_data(key, value):
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


def get_data_pyodide(key, globals_dict=None):
    """
    Request data from the host environment through a JavaScript function defined in the
    JupyterLite extension `data_bridge`.
    Args:
        key (string): The name under which data is expected to be received.
        globals_dict (dict): A dictionary to store the received data. Defaults to None.
    """
    js_code = f"""
    (function() {{
        if (window.requestDataFromHost) {{
            window.requestDataFromHost('{key}')
        }} else {{
            console.error('requestDataFromHost function is not defined on the window object.')
        }}
    }})();
    """
    display(Javascript(js_code))
    print(f"Status: {key} requested")


def get_data_python(key, globals_dict=None):
    """
    Read data from the `uploads` folder in a JupyterLab environment.
    Args:
        key (string): The name under which data is expected to be received.
        globals_dict (dict): A dictionary to store the received data. Defaults to None.
    """
    try:
        materials = []
        for filename in os.listdir(UPLOADS_FOLDER):
            if filename.endswith(".json"):
                with open(os.path.join(UPLOADS_FOLDER, filename), "r") as file:
                    data = json.load(file)
                name = os.path.splitext(filename)[0]
                print(f"Data from {name} has been read successfully.")
                materials.append(data)
        if globals_dict is not None:
            globals_dict[key] = materials
    except FileNotFoundError:
        print("No data found in the 'uploads' folder.")


def get_data(key, globals_dict=None):
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
