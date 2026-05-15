import re
from typing import List

from src.py.mat3ra.notebooks_utils.primitive.environment import ENVIRONMENT


def install_package_python(pkg: str, verbose: bool = True):
    """
    Install a package in a standard Python environment.

    Args:
        pkg (str): The name of the package to install.
        verbose (bool): Whether to print the name of the installed package.
    """
    # NOTE: in a regular Python environment packages should be installed via pip,
    # not programmatically from config.yml. Direct user to do so.
    print('To install packages, run `pip install ".[all]"` in the terminal')


def get_packages_list(requirements_dict: dict, notebook_name_pattern: str = "") -> List[str]:
    """
    Get the list of packages to install based on the requirements dict.

    Args:
        requirements_dict (dict): The dictionary containing the requirements.
        notebook_name_pattern (str): The pattern of the notebook name.

    Returns:
        List[str]: The list of packages to install.
    """
    packages_default_common = requirements_dict.get("default", {}).get("packages_common", [])
    packages_default_environment_specific = requirements_dict.get("default", {}).get(
        f"packages_{ENVIRONMENT.value}", []
    )

    matching_notebook_requirements_list = [
        cfg for cfg in requirements_dict.get("notebooks", []) if re.search(cfg.get("name"), notebook_name_pattern)
    ]
    packages_notebook_common = []
    packages_notebook_environment_specific = []

    for notebook_requirements in matching_notebook_requirements_list:
        packages_common = notebook_requirements.get("packages_common", [])
        packages_environment_specific = notebook_requirements.get(f"packages_{ENVIRONMENT.value}", [])
        if packages_common:
            packages_notebook_common.extend(packages_common)
        if packages_environment_specific:
            packages_notebook_environment_specific.extend(packages_environment_specific)

    # Note: environment specific packages have to be installed first,
    # because in Pyodide common packages might depend on them
    return [
        *packages_default_environment_specific,
        *packages_notebook_environment_specific,
        *packages_default_common,
        *packages_notebook_common,
    ]
