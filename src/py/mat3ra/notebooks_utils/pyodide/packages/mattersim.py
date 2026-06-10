import sys
from importlib import import_module

import numpy as np
import torch

from ...primitive.environment import is_pyodide_environment
from .torch import _make_stub_module, patch_torch_distributed


def apply_patches():
    patch_torch_distributed()
    patch_mattersim_deps()


def patch_mattersim_deps():
    """
    Stub heavy dependencies required by MatterSim but not needed for inference.

    Stubs: loguru, azure.*, atomate2, seekpath, phonopy, phono3py, mp_api,
    sklearn, and patches e3nn to disable JIT/torch.compile.
    """
    loguru_mod = _make_stub_module("loguru")

    class _Logger:
        def info(self, msg, *a, **k):
            print(f"INFO: {msg}")

        def warning(self, msg, *a, **k):
            print(f"WARNING: {msg}")

        def error(self, msg, *a, **k):
            print(f"ERROR: {msg}")

        def debug(self, msg, *a, **k):
            pass

        def trace(self, msg, *a, **k):
            pass

        def success(self, msg, *a, **k):
            print(f"✓ {msg}")

        def __getattr__(self, name):
            return lambda *a, **k: None

    loguru_mod.logger = _Logger()

    _make_stub_module("azure", submodules=["identity", "storage", "storage.blob"])

    for package_name in [
        "atomate2",
        "seekpath",
        "phonopy",
        "phono3py",
        "mp_api",
        "jobflow",
        "emmet",
        "emmet.core",
        "emmet.core.tasks",
        "maggma",
    ]:
        _make_stub_module(package_name)

    _patch_sklearn()
    _patch_pyodide_http()
    _patch_e3nn()
    _patch_torch_jit()
    _patch_torch_ema()
    _patch_torchmetrics()
    _patch_torch_geometric()
    _patch_torch_runstats()

    print("✓ MatterSim dependency stubs applied")


def _patch_sklearn():
    _make_stub_module(
        "sklearn",
        submodules=[
            "base",
            "utils",
            "utils.validation",
            "preprocessing",
            "model_selection",
            "gaussian_process",
            "gaussian_process.kernels",
        ],
    )

    class _GPR:
        def __init__(self, *a, **k):
            pass

        def fit(self, *a, **k):
            return self

        def predict(self, X, return_std=False):
            mean = np.zeros(X.shape[0])
            if return_std:
                return mean, np.ones(X.shape[0])
            return mean

        def log_marginal_likelihood(self):
            return 0.0

    sys.modules["sklearn.gaussian_process"].GaussianProcessRegressor = _GPR

    class _Kernel:
        pass

    class _DotProduct(_Kernel):
        def __init__(self, *a, **k):
            pass

    class _Hyperparameter:
        def __init__(self, *a, **k):
            pass

    sklearn_kernels = sys.modules["sklearn.gaussian_process.kernels"]
    sklearn_kernels.Kernel = _Kernel
    sklearn_kernels.DotProduct = _DotProduct
    sklearn_kernels.Hyperparameter = _Hyperparameter


def _patch_pyodide_http():
    try:
        import pyodide_http

        pyodide_http.patch_all()
    except ImportError:
        pass


def _patch_e3nn():
    try:
        import e3nn

        e3nn._SO3_INITIALIZED = True
    except Exception:
        pass


def _patch_torch_jit():
    if not hasattr(torch.jit, "_original_script"):

        def _noop_script(obj=None, *a, **k):
            if obj is not None:
                return obj
            return lambda fn: fn

        torch.jit.script = _noop_script


def _patch_torch_ema():
    torch_ema_mod = _make_stub_module("torch_ema")

    class _EMA:
        def __init__(self, *a, **k):
            pass

    torch_ema_mod.ExponentialMovingAverage = _EMA


def _patch_torchmetrics():
    torchmetrics_mod = _make_stub_module("torchmetrics")

    class _MeanMetric:
        def __init__(self, *a, **k):
            pass

    torchmetrics_mod.MeanMetric = _MeanMetric


def _patch_torch_geometric():
    _make_stub_module("torch_geometric", submodules=["data", "loader", "utils"])

    class _Data:
        def __init__(self, **kwargs):
            for k, v in kwargs.items():
                setattr(self, k, v)

        def to(self, device):
            for k, v in self.__dict__.items():
                if isinstance(v, torch.Tensor):
                    setattr(self, k, v.to(device))
            return self

    sys.modules["torch_geometric.data"].Data = _Data

    class _DataLoader:
        def __init__(self, dataset, batch_size=1, shuffle=False, **kwargs):
            self._dataset = list(dataset)

        def __iter__(self):
            for item in self._dataset:
                for attr in list(vars(item).keys()):
                    val = getattr(item, attr)
                    if isinstance(val, (int, float)):
                        setattr(item, attr, torch.tensor([val]))
                if not hasattr(item, "num_graphs"):
                    item.num_graphs = 1
                if not hasattr(item, "batch"):
                    n_atoms = item.num_atoms if hasattr(item, "num_atoms") else torch.tensor([0])
                    if isinstance(n_atoms, torch.Tensor):
                        n_atoms = int(n_atoms.item())
                    item.batch = torch.zeros(n_atoms, dtype=torch.long)
                yield item

        def __len__(self):
            return len(self._dataset)

    sys.modules["torch_geometric.loader"].DataLoader = _DataLoader


def _patch_torch_runstats():
    _make_stub_module("torch_runstats", submodules=["scatter"])

    def _scatter(src, index, dim_size=None, dim=0, reduce="sum"):
        if dim_size is None:
            dim_size = int(index.max()) + 1
        out = torch.zeros(dim_size, *src.shape[1:], dtype=src.dtype, device=src.device)
        if src.dim() == 1:
            idx = index
        else:
            idx = index.unsqueeze(-1).expand_as(src)
        if reduce == "sum" or reduce == "add":
            out.scatter_add_(0, idx, src)
        elif reduce == "mean":
            out.scatter_add_(0, idx, src)
            count = torch.zeros(dim_size, dtype=src.dtype, device=src.device)
            count.scatter_add_(0, index, torch.ones(index.shape[0], dtype=src.dtype, device=src.device))
            count = count.clamp(min=1)
            if src.dim() > 1:
                count = count.unsqueeze(-1)
            out = out / count
        return out

    sys.modules["torch_runstats.scatter"].scatter = _scatter

    def _scatter_mean(src, index, dim_size=None, dim=0):
        return _scatter(src, index, dim_size=dim_size, dim=dim, reduce="mean")

    sys.modules["torch_runstats.scatter"].scatter_mean = _scatter_mean


MODEL_PATHS_MAP = {
    "1m": "/drive/packages/models/mattersim-v1.0.0-1M.pth",
}


def get_mattersim_model_pyodide(model: str, device="cpu", **kwargs):
    if model not in MODEL_PATHS_MAP:
        raise ValueError(f"Invalid model name: {model}. Valid options are: {list(MODEL_PATHS_MAP.keys())}")
    forcefield = import_module("mattersim.forcefield")
    return forcefield.MatterSimCalculator.from_checkpoint(load_path=MODEL_PATHS_MAP[model], device=device, **kwargs)


def create_mattersim_calculator(model="1m", device="cpu", model_path=None, checkpoint=None, **kwargs):
    if is_pyodide_environment():
        return get_mattersim_model_pyodide(model=model, device=device, **kwargs)

    resolved_model_path = model_path or checkpoint

    forcefield = import_module("mattersim.forcefield")
    return forcefield.MatterSimCalculator.from_checkpoint(load_path=str(resolved_model_path), device=device, **kwargs)
