from importlib import import_module

from .torch import apply_common_torch_patches

MLFF_PATCH_MODULES = {
    "mace": "mat3ra.notebooks_utils.pyodide.packages.mace",
    "uma": "mat3ra.notebooks_utils.pyodide.packages.uma",
    "mattersim": "mat3ra.notebooks_utils.pyodide.packages.mattersim",
    "nequip": "mat3ra.notebooks_utils.pyodide.packages.nequip",
}


def normalize_mlff_name(mlff_name: str) -> str:
    return (mlff_name or "").strip().lower()


def apply_mlff_patches(mlff_name: str):
    mlff = normalize_mlff_name(mlff_name)
    if mlff not in MLFF_PATCH_MODULES:
        raise ValueError(f"Unsupported MLFF: {mlff_name!r}")

    patch_module = import_module(MLFF_PATCH_MODULES[mlff])
    patch_module.apply_patches()


def apply_all_patches(mlff_name: str):
    apply_common_torch_patches()
    apply_mlff_patches(mlff_name)
    print("\n✅ All Pyodide patches applied successfully!")
