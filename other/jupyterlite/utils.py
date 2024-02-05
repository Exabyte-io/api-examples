from IPython.display import display, Javascript
import json
import time

"""
This module contains a function to install packages in a Pyodide environment.
Pyodide uses micropip as replacement for pip to install packages.
Package must be compiled for none-any platform.
"""

try:
    import micropip
except ImportError:
    raise ImportError(
        "This module intended to be used in a Pyodide environment. Please install packages ypurself using pip."
    )


async def install_packages(notebook_name, requirements_path="config.yml", verbose=True):
    await micropip.install("pyyaml")
    import yaml

    with open(requirements_path, "r") as f:
        requirements = yaml.safe_load(f)

    packages = None
    for requirement in requirements:
        if requirement["notebook"] == notebook_name:
            packages = requirement["packages"]
            break

    if packages is None:
        raise ValueError(f"No packages found for notebook {notebook_name}")

    async def install_package(pkg):
        """
        Installs a package in a Pyodide environment.
        Args:
            pkg: The name of the package to install.

        Returns:
            None
        """
        is_url = pkg.startswith("http://") or pkg.startswith("https://")
        are_dependencies_installed = not is_url
        await micropip.install(pkg, deps=are_dependencies_installed)
        # Extract package name for printing
        pkg_name = pkg.split("/")[-1].split("-")[0] if is_url else pkg.split("==")[0]
        if verbose:
            print(f"Installed {pkg_name}")

    for package in packages:
        await install_package(package)
    if verbose:
        print("All packages installed.")


def set_data(key, value):
    """
    This function takes a Python object, serializes it to JSON, and sends it to the host environment
    through a JavaScript function defined in the JupyterLite extension `data_bridge`.

    Args:
        materials (object): The Python object to be sent to the host environment.
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
    This function requests materials from the host environment through a JavaScript function defined in the JupyterLite
    extension `data_bridge`. The materials are then returned to the Python environment.
    """
    js_code = """
    (function() {
        if (window.requestDataFromHost) {
            window.requestDataFromHost();
            
        } else {
            console.error('requestDataFromHost function is not defined on the window object.');
        }
    })();
    """

    display(Javascript(js_code))
    time.sleep(3)  # JS postMessage is asynchronous, so we need to wait for the response from JS host
    print(f"Status: {key} received")
