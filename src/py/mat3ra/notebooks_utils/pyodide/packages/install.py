import importlib
import os
import sys

from ...primitive.logger import log

try:
    import micropip  # type: ignore
except ImportError:
    micropip = None  # type: ignore

PYODIDE_INIT_PACKAGES = ["pyyaml"]
PYODIDE_INIT_MODULES = ["yaml"]


def get_config_yml_file_path(config_file_path: str) -> str:
    """
    Resolve the absolute path to the config.yml file.

    Args:
        config_file_path (str): Relative path override; empty string uses the JupyterLite default (/drive/config.yml).

    Returns:
        str: Absolute path to the config file.
    """
    config_file_full_path = os.path.normpath(os.path.join("/drive/", "./config.yml"))
    if config_file_path != "":
        config_file_full_path = os.path.normpath(os.path.join(os.getcwd(), config_file_path))
    return config_file_full_path


async def install_init():
    if sys.platform != "emscripten":
        return

    await micropip.install(PYODIDE_INIT_PACKAGES)
    for module in PYODIDE_INIT_MODULES:
        importlib.import_module(module)


async def install_package_pyodide(pkg: str, verbose: bool = True):
    """
    Install a package in a Pyodide environment.

    Args:
        pkg (str): The name of the package to install. Can be prefixed with 'nodeps:' to skip dependencies.
        verbose (bool): Whether to print the name of the installed package.

    Examples:
        await install_package_pyodide("numpy")  # installs with deps
        await install_package_pyodide("nodeps:e3nn==0.4.4")  # installs without deps
    """
    if pkg.startswith("nodeps:"):
        pkg = pkg.replace("nodeps:", "")
        are_dependencies_installed = False
    else:
        is_url = pkg.startswith("http://") or pkg.startswith("https://") or pkg.startswith("emfs:/")
        are_dependencies_installed = not is_url

    await micropip.install(pkg, deps=are_dependencies_installed)
    pkg_name = pkg.split("/")[-1].split("-")[0] if "://" in pkg else pkg.split("==")[0]
    if verbose:
        log(f"Installed {pkg_name}", force_verbose=verbose)
