from IPython.display import display, Javascript
import json
import os
from enum import Enum

UPLOADS_FOLDER = "uploads"


class Environment(Enum):
    PYODIDE = "pyodide"
    PYTHON = "python"


# Environment detection
# default value for env.HOME from https://pyodide.org/en/stable/usage/api/js-api.html
ENVIRONMENT = Environment.PYODIDE if os.environ.get("HOME") == "/home/pyodide" else Environment.PYTHON

if ENVIRONMENT == Environment.PYODIDE:
    import micropip

if ENVIRONMENT == Environment.PYTHON:
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
    subprocess.check_call([sys.executable, "-m", "pip", "install", "--upgrade", "pip"])
    subprocess.check_call([sys.executable, "-m", "pip", "install", pkg])
    if verbose:
        print(f"Installed {pkg}")


async def install_packages(notebook_name, requirements_path="config.yml", verbose=True):
    """
    Install the packages listed in the requirements file for the notebook with the given name.
    Args:
        notebook_name (string): The name of the notebook for which to install packages.
        requirements_path (string): The path to the requirements file.
        verbose (bool): Whether to print the names of the installed packages and status of installation.
    """
    if ENVIRONMENT == Environment.PYODIDE:
        await micropip.install("pyyaml")
    if ENVIRONMENT == Environment.PYTHON:
        install_package_python("pyyaml")
    import yaml

    with open(requirements_path, "r") as f:
        requirements = yaml.safe_load(f)

    default_packages = requirements.get("default", {})
    notebook_packages_common = None
    notebook_packages = None
    for notebook in requirements.get("notebooks", []):
        if notebook.get("name") == notebook_name:
            notebook_packages_common = notebook.get("packages", [])
            packages_key = "packages_pyodide" if ENVIRONMENT == Environment.PYODIDE else "packages_python"
            notebook_packages = notebook.get(packages_key, [])
            break

    # Hash the requirements to avoid re-installing packages
    requirements_hash = str(hash(json.dumps(requirements)))
    if os.environ.get("requirements_hash") != requirements_hash:
        if verbose:
            print("Installing packages...")
        for pkg in default_packages:
            if ENVIRONMENT == Environment.PYODIDE:
                await install_package_pyodide(pkg, verbose)
            if ENVIRONMENT == Environment.PYTHON:
                install_package_python(pkg, verbose)
        if notebook_packages_common:
            for pkg in notebook_packages_common:
                if ENVIRONMENT == Environment.PYODIDE:
                    await install_package_pyodide(pkg, verbose)
                if ENVIRONMENT == Environment.PYTHON:
                    install_package_python(pkg, verbose)
        if notebook_packages:
            for pkg in notebook_packages:
                if ENVIRONMENT == Environment.PYODIDE:
                    await install_package_pyodide(pkg, verbose)
                if ENVIRONMENT == Environment.PYTHON:
                    install_package_python(pkg, verbose)
        if verbose:
            print("Packages installed.")

    os.environ["requirements_hash"] = requirements_hash


def set_data(key, value):
    """
    This function takes a Python object, serializes it to JSON, and sends it to the host environment
    through a JavaScript function defined in the JupyterLite extension `data_bridge`.

    Args:
        key (string): The name under which data will be sent.
        value (Any): The value to send to the host environment.
    """
    if ENVIRONMENT == Environment.PYODIDE:
        # Pyodide: Send data to host environment using JavaScript
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
    elif ENVIRONMENT == Environment.PYTHON:
        # Standard Python environment: Write data to 'uploads' folder
        if not os.path.exists(UPLOADS_FOLDER):
            os.makedirs(UPLOADS_FOLDER)
        file_path = os.path.join(UPLOADS_FOLDER, f"{key}.json")
        with open(file_path, "w") as file:
            json.dump(value, file)
        print(f"Data for {key} written to {file_path}")


def get_data(key, globals_dict=None):
    """
    Request data from the host environment through a JavaScript function defined in the
    JupyterLite extension `data_bridge` or read the data directly from the `uploads` folder in a JupyterLab environment.
    Args:
        key (string): The name under which data is expected to be received.
        globals_dict (dict): A dictionary to store the received data. Defaults to None.
    """
    if ENVIRONMENT == Environment.PYODIDE:
        # JupyterLite environment: Request data using JavaScript extension
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
    elif ENVIRONMENT == Environment.PYTHON:
        # JupyterLab environment: Read data from the 'uploads' folder
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
