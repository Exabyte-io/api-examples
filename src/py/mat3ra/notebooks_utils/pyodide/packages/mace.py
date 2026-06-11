import sys
import types
from importlib import import_module

from ...primitive.environment import is_pyodide_environment


def apply_patches():
    patch_mace_training()
    patch_mace_tools()


def patch_mace_training():
    """
    Stub lmdb and h5py packages.

    These are C-extension packages used by MACE's training/dataset code
    but not needed for inference. Stubs allow imports to succeed.
    """
    for package_name in ("lmdb", "h5py"):
        if package_name not in sys.modules:
            sys.modules[package_name] = types.ModuleType(package_name)

    print("✓ LMDB and HDF5 stubs applied")


def patch_mace_tools():
    """
    Fix MACE's torch_geometric import order issues in Pyodide.

    In Pyodide, torch_geometric.data may not be set during circular imports.
    Pre-importing ensures the attribute is available when MACE needs it.
    """
    try:
        torch_geometric = import_module("mace.tools.torch_geometric")
        torch_geometric_data = import_module("mace.tools.torch_geometric.data")
        torch_geometric.data = torch_geometric_data
        print("✓ MACE tools patches applied")
    except Exception as exc:
        print(f"⚠ MACE tools patches skipped: {exc}")


MODEL_PATHS_MAP = {
    "small": "/drive/packages/models/2023-12-10-mace-128-L0_energy_epoch-249.model",
    "medium": "/drive/packages/models/2023-12-03-mace-128-L1_epoch-199.model",
    "large": "/drive/packages/models/MACE_MPtrj_2022.9.model",
}


def get_mace_model_pyodide(model: str, dispersion=False, default_dtype="float32", device="cpu", **kwargs):
    if model not in MODEL_PATHS_MAP:
        raise ValueError(f"Invalid model name: {model}. Valid options are: {list(MODEL_PATHS_MAP.keys())}")
    model_path = MODEL_PATHS_MAP[model]
    mace_calculators = import_module("mace.calculators")
    return mace_calculators.MACECalculator(
        model_path=model_path, dispersion=dispersion, default_dtype=default_dtype, device=device, **kwargs
    )


def create_mace_calculator(model="large", dispersion=True, default_dtype="float32", device="cpu", **kwargs):
    if is_pyodide_environment():
        return get_mace_model_pyodide(
            model=model,
            dispersion=dispersion,
            default_dtype=default_dtype,
            device=device,
            **kwargs,
        )

    mace_calculators = import_module("mace.calculators")
    return mace_calculators.mace_mp(
        model=model,
        dispersion=dispersion,
        default_dtype=default_dtype,
        device=device,
        **kwargs,
    )
