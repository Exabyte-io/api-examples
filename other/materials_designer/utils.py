"""
This module contains a function to install packages in a Pyodide environment.
Pyodide uses micropip as replacement for pip to install packages.
Package must be compiled for none-any platform.
"""
import micropip
import json


async def install_packages(verbose=True):
    packages = json.loads(open("packages.json").read())["packages"]

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
