from IPython.display import display, Javascript
import json
import os

try:
    import micropip
except ImportError:
    raise ImportError(
        "This module intended to be used in a Pyodide environment. Please install packages yourself using pip."
    )

async def install_package(pkg, verbose=True):
        """
        Installs a package in a Pyodide environment.
        Args:
            pkg: The name of the package to install.
            verbose: Whether to print the name of the installed package.

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

async def install_packages(notebook_name, requirements_path="config.yml", verbose=True):
    """
    This function installs the packages listed in the requirements file for the notebook with the given name.
    Args:
        notebook_name: The name of the notebook for which to install packages.
        requirements_path: The path to the requirements file.
        verbose: Whether to print the names of the installed packages and status of installation.
    """
    await micropip.install("pyyaml")
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
            await install_package(package)
        if verbose:
            print("All packages installed.")
    
    os.environ["requirements_hash"] = requirements_hash


def set_data(key, value):
    """
    This function takes a Python object, serializes it to JSON, and sends it to the host environment
    through a JavaScript function defined in the JupyterLite extension `data_bridge`.

    Args:
        key: The key to use for the data.
        value: The value to send to the host environment.
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
    This function requests data from the host environment through a JavaScript function defined in the JupyterLite
    extension `data_bridge`. The data is then returned to the Python environment.
    Args:
        key: The key to use for the data.
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
