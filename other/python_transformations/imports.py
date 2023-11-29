"""Installation of all supported packages"""
import micropip

for pkg in [
    "https://files.mat3ra.com:44318/uploads/pymatgen-2023.9.10-py3-none-any.whl",
    "https://files.mat3ra.com:44318/web/pyodide/spglib-2.0.2-py3-none-any.whl",
    "https://files.pythonhosted.org/packages/d9/0e/2a05efa11ea33513fbdf4a2e2576fe94fd8fa5ad226dbb9c660886390974/ruamel.yaml-0.17.32-py3-none-any.whl",
]:
    await micropip.install(pkg, deps=False)
    print("Installed " + pkg)

for pkg in [
    "ase==3.22.1",
    "networkx==3.1",
    "monty==2023.11.3",
    "scipy==1.11.1",
    "lzma",
    "tabulate==0.9.0",
    "sqlite3",
    "sympy==1.12",
]:
    await micropip.install(pkg)
    print("Installed " + pkg)
