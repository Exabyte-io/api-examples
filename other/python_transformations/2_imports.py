"""TITLE: Import of all supported packages."""
"""BLOCK: Install Packages"""
"""
This scripts installs the required packages for a basic Pyodide environment for Materials structure manipulation.
It uses micropip - a Python package manager for Pyodide.
"""
import micropip

packages = [
    "https://files.mat3ra.com:44318/uploads/pymatgen-2023.9.10-py3-none-any.whl",
    "https://files.mat3ra.com:44318/web/pyodide/spglib-2.0.2-py3-none-any.whl",
    "https://files.pythonhosted.org/packages/d9/0e/2a05efa11ea33513fbdf4a2e2576fe94fd8fa5ad226dbb9c660886390974/ruamel.yaml-0.17.32-py3-none-any.whl",
    "ase==3.22.1",
    "networkx==3.2.1",
    "monty==2023.11.3",
    "scipy==1.11.1",
    "lzma",
    "tabulate==0.9.0",
    "sqlite3",
    "sympy==1.12",
]


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
    print(f"Installed {pkg_name}")


for package in packages:
    await install_package(package)
