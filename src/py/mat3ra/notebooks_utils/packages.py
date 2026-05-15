import importlib
import json
import os
from typing import List

from src.py.mat3ra.notebooks_utils.ipython.packages.install import get_packages_list, install_package_python

from .primitive.enums import EnvironmentsEnum
from .primitive.environment import ENVIRONMENT, is_pyodide_environment
from .primitive.logger import log
from .pyodide.packages.install import get_config_yml_file_path, install_init


async def install_package(pkg: str, verbose: bool = True):
    """
    Install a package in the current environment.

    Args:
        pkg (str): The name of the package to install.
        verbose (bool): Whether to print the name of the installed package.
    """
    if ENVIRONMENT == EnvironmentsEnum.PYODIDE:
        from .pyodide.packages.install import install_package_pyodide

        await install_package_pyodide(pkg, verbose)
    elif ENVIRONMENT == EnvironmentsEnum.PYTHON:
        install_package_python(pkg, verbose)


async def install_packages_with_hashing(packages: List[str], verbose: bool = True):
    """
    Install the given packages, skipping if the list is unchanged since last run.

    Args:
        packages (List[str]): The list of packages to install.
        verbose (bool): Whether to print the names of the installed packages.
    """
    requirements_hash = str(hash(json.dumps(packages)))
    if os.environ.get("requirements_hash") != requirements_hash:
        for pkg in packages:
            await install_package(pkg, verbose)
        if verbose:
            log("Packages installed successfully.", force_verbose=verbose)
        os.environ["requirements_hash"] = requirements_hash
    else:
        if verbose:
            log("Packages are already installed.", force_verbose=verbose)


async def install_packages(notebook_name_pattern: str, config_file_path: str = "", verbose: bool = True):
    """
    Install the packages listed in config.yml for the given notebook name pattern.

    Usage in notebooks:
        from mat3ra.notebooks_utils.packages import install_packages
        await install_packages("my_notebook")

    Args:
        notebook_name_pattern (str): Pattern matched against notebook names in config.yml.
        config_file_path (str): Path to config.yml; empty string uses the JupyterLite default (/drive/config.yml).
        verbose (bool): Whether to print install progress.
    """
    if is_pyodide_environment():
        await install_init()

    yaml = importlib.import_module("yaml")
    with open(get_config_yml_file_path(config_file_path), "r") as f:
        requirements_dict = yaml.safe_load(f)

    packages = get_packages_list(requirements_dict, notebook_name_pattern)
    await install_packages_with_hashing(packages, verbose)
