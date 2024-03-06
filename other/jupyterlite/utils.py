from IPython.display import display, Javascript
import json
import os

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
    Installs a package in a Pyodide environment.
    """
    is_url = pkg.startswith("http://") or pkg.startswith("https://")
    are_dependencies_installed = not is_url
    await micropip.install(pkg, deps=are_dependencies_installed)
    pkg_name = pkg.split("/")[-1].split("-")[0] if is_url else pkg.split("==")[0]
    if verbose:
        print(f"Installed {pkg_name}")


def install_package_python(pkg, verbose=True):
    """
    Installs a package in a standard Python environment.
    """
    subprocess.check_call([sys.executable, "-m", "pip", "install", pkg])
    if verbose:
        print(f"Installed {pkg}")


async def install_packages(notebook_name, requirements_path="config.yml", verbose=True):
    """
    Installs the packages listed in the requirements file for the notebook with the given name.
    """
    if IN_PYODIDE:
        await micropip.install("pyyaml")
    else:
        install_package_python("pyyaml", verbose=verbose)

    import yaml

    requirements_hash = ""

    with open(requirements_path, "r") as f:
        requirements = yaml.safe_load(f)
        requirements_hash = str(hash(json.dumps(requirements)))

    packages = None
    for requirement in requirements:
        if requirement["notebook"] == notebook_name:
            packages = requirement["packages"]
            break

    if packages is None:
        raise ValueError(f"No packages found for notebook {notebook_name}")

    if os.environ.get("requirements_hash") != requirements_hash:
        for package in packages:
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
    serialized_data = json.dumps({key: value})
    js_code = f"""
    (function() {{
        window.sendDataToHost({serialized_data})
        console.log({serialized_data})
    }})();
    """

    display(Javascript(js_code))
    print(f"Status: {key} sent to host.")


def get_data(key):
    """
    This function either requests data from the host environment through a JavaScript function defined in the
    JupyterLite extension `data_bridge` or reads the data directly from the `uploads` folder in a JupyterLab environment.
    Args:
        key (string): The name under which data is expected to be received or the file name to read in JupyterLab.
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
            uploads_path = 'uploads'
            materials = []
            for filename in os.listdir(uploads_path):
                if filename.endswith('.json'):
                    with open(os.path.join(uploads_path, filename), 'r') as file:
                        data = json.load(file)
                    key = os.path.splitext(filename)[0]
                    print(f"Data from {key} has been read successfully.")
                    materials.append(data)
            globals()[key] = materials

        except FileNotFoundError:
            print(f"File {key} not found in 'uploads' folder.")

