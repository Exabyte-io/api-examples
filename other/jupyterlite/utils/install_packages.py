"""
This module contains a function to install packages in a Pyodide environment.
Pyodide uses micropip as replacement for pip to install packages.
Package must be compiled for none-any platform.
"""

import yaml

try:
    import micropip
except ImportError:
    raise ImportError(
        "This module intended to be used in a Pyodide environment. Please install packages ypurself using pip."
    )


async def install_packages(notebook_name, requirements_path="requirements.yml", verbose=True):
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
