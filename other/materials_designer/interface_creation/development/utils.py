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

import json

packages = [
    "https://files.mat3ra.com:44318/uploads/pymatgen-2023.9.10-py3-none-any.whl",
    "https://files.mat3ra.com:44318/web/pyodide/spglib-2.0.2-py3-none-any.whl",
    "https://files.pythonhosted.org/packages/d9/0e/2a05efa11ea33513fbdf4a2e2576fe94fd8fa5ad226dbb9c660886390974/ruamel.yaml-0.17.32-py3-none-any.whl",
    "ase==3.22.1",
    "networkx==3.2.1",
    "monty==2023.11.3",
    "scipy==1.11.2",
    "lzma",
    "tabulate==0.9.0",
    "sqlite3",
    "sympy==1.12",
    "uncertainties==3.1.6",
    "ipywidgets",
]


async def install_packages(verbose=True):

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
        print("Done!")