"""Installation of all supported packages"""
import micropip

await micropip.install("https://files.mat3ra.com:44318/uploads/pymatgen-2023.9.10-py3-none-any.whl", deps=False)
await micropip.install("https://files.mat3ra.com:44318/web/pyodide/spglib-2.0.2-py3-none-any.whl", deps=False)
await micropip.install(
    "https://files.pythonhosted.org/packages/d9/0e/2a05efa11ea33513fbdf4a2e2576fe94fd8fa5ad226dbb9c660886390974/ruamel.yaml-0.17.32-py3-none-any.whl",
    deps=False,
)
for pkg in ["ase", "networkx", "monty", "scipy", "lzma", "tabulate", "sqlite3", "sympy"]:
    await micropip.install(pkg)
