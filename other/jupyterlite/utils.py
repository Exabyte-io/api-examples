from IPython.display import display, Javascript
import json
import os

UPLOADS_FOLDER = "uploads"

# Environment detection
IN_PYODIDE = False

try:
    import micropip

    IN_PYODIDE = True
except ImportError:
    IN_PYODIDE = False

if not IN_PYODIDE:
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
    if IN_PYODIDE:
        await micropip.install("pyyaml")
    else:
        install_package_python("pyyaml", verbose=verbose)

    import yaml

    with open(requirements_path, "r") as f:
        requirements = yaml.safe_load(f)
    # Hash the requirements to avoid re-installing packages
    requirements_hash = str(hash(json.dumps(requirements)))

    default_packages = requirements.get("default", {}).get("packages", [])

    # Install packages in Pyodide that are loaded by default in Python
    for package in default_packages:
        if IN_PYODIDE:
            await install_package_pyodide(package, verbose=verbose)

    notebook_packages = None
    for notebook in requirements.get("notebooks", []):
        if notebook.get("notebook") == notebook_name:
            notebook_packages = notebook.get("packages", [])
            break

    if notebook_packages is None:
        raise ValueError(f"No packages found for notebook {notebook_name}")

    if os.environ.get("requirements_hash") != requirements_hash:
        for package in notebook_packages:
            if IN_PYODIDE:
                await install_package_pyodide(package, verbose=verbose)
            else:
                install_package_python(package, verbose=verbose)
        if verbose:
            print("All packages installed.")

    os.environ["requirements_hash"] = requirements_hash


def set_data(key, value):
    """
    This function takes a Python object, serializes it to JSON, and sends it to the host environment
    through a JavaScript function defined in the JupyterLite extension `data_bridge`.

    Args:
        key (string): The name under which data will be sent.
        value (Any): The value to send to the host environment.
    """
    if IN_PYODIDE:
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
    else:
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
    if IN_PYODIDE:
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
    else:
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
