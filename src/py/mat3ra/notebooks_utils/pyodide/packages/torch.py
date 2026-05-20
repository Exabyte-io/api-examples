"""
Patches for PyTorch and related packages to work in Pyodide environment.

This module provides patches for various torch-related packages that don't work
in Pyodide's WASM environment, organized by functionality.

Usage:
    from mat3ra.notebooks_utils.other.torch_pyodide import (
        patch_torch_linalg,
        patch_torch_testing,
        patch_matscipy,
        patch_lmdb_h5py,
        patch_mace_tools,
    )

    # Apply all patches
    patch_torch_linalg()
    patch_torch_testing()
    patch_matscipy()
    patch_lmdb_h5py()
    patch_mace_tools()
"""

import sys
import types
from collections import namedtuple

import numpy as np
import torch

# Define return types to mimic PyTorch's named tuples
EigRet = namedtuple("linalg_eig", ["eigenvalues", "eigenvectors"])  # type: ignore
EighRet = namedtuple("linalg_eigh", ["eigenvalues", "eigenvectors"])  # type: ignore
LUFactorReturn = namedtuple("LUFactorReturn", ["LU", "pivots"])


def _to_np(tensor):
    return tensor.detach().cpu().numpy()


def _to_torch(array, device, dtype=None):
    if dtype:
        return torch.tensor(array, device=device, dtype=dtype)
    return torch.tensor(array, device=device)


# ==============================================================================
# Torch linalg patches
# ==============================================================================


def _patch_solve(A, B, *args, **kwargs):
    return _to_torch(np.linalg.solve(_to_np(A), _to_np(B)), A.device)


def _patch_inv(A, *args, **kwargs):
    return _to_torch(np.linalg.inv(_to_np(A)), A.device)


def _patch_det(A, *args, **kwargs):
    return _to_torch(np.linalg.det(_to_np(A)), A.device)


def _patch_cholesky(A, *args, **kwargs):
    return _to_torch(np.linalg.cholesky(_to_np(A)), A.device)


def _patch_eig(A, *args, **kwargs):
    vals, vecs = np.linalg.eig(_to_np(A))
    return EigRet(_to_torch(vals, A.device), _to_torch(vecs, A.device))


def _patch_eigh(A, UPLO="L", *args, **kwargs):
    vals, vecs = np.linalg.eigh(_to_np(A), UPLO=UPLO)
    return EighRet(_to_torch(vals, A.device), _to_torch(vecs, A.device))


def _patch_lu_factor(A, *args, **kwargs):
    import scipy.linalg

    A_np = _to_np(A)
    if A_np.ndim > 2:
        orig_shape = A_np.shape
        A_reshaped = A_np.reshape(-1, orig_shape[-2], orig_shape[-1])
        lu_list, piv_list = [], []
        for mat in A_reshaped:
            lu, piv = scipy.linalg.lu_factor(mat)
            lu_list.append(lu)
            piv_list.append(piv)
        LU_np = np.stack(lu_list).reshape(orig_shape)
        piv_np = np.stack(piv_list).reshape(orig_shape[:-2] + (-1,))
    else:
        LU_np, piv_np = scipy.linalg.lu_factor(A_np)
    return LUFactorReturn(_to_torch(LU_np, A.device, A.dtype), _to_torch(piv_np, A.device, torch.int32))


def _patch_lu_solve(LU, pivots, B, *args, **kwargs):
    import scipy.linalg

    LU_np, piv_np, B_np = _to_np(LU), _to_np(pivots), _to_np(B)
    if LU_np.ndim > 2:
        orig_shape_B = B_np.shape
        LU_reshaped = LU_np.reshape(-1, LU_np.shape[-2], LU_np.shape[-1])
        piv_reshaped = piv_np.reshape(-1, piv_np.shape[-1])
        is_vector = B_np.ndim == LU_np.ndim - 1
        B_reshaped = (
            B_np.reshape(-1, B_np.shape[-1], 1) if is_vector else B_np.reshape(-1, B_np.shape[-2], B_np.shape[-1])
        )
        X_list = [scipy.linalg.lu_solve((lu, piv), b) for lu, piv, b in zip(LU_reshaped, piv_reshaped, B_reshaped)]
        X_np = np.stack(X_list).reshape(orig_shape_B)
    else:
        X_np = scipy.linalg.lu_solve((LU_np, piv_np), B_np)
    return _to_torch(X_np, B.device, B.dtype)


def _tensor_array_compat(self, dtype=None):
    """Replacement for Tensor.__array__ in Pyodide where tensor.numpy() is unavailable."""
    arr = np.array(self.tolist())
    return arr.astype(dtype) if dtype is not None else arr


def patch_torch_linalg():
    """
    Patch torch.linalg functions to use NumPy/SciPy implementations.

    This fixes LAPACK functions that are not available in Pyodide's WASM torch build.
    """
    torch.linalg.solve = _patch_solve
    torch.linalg.inv = _patch_inv
    torch.inverse = _patch_inv  # Alias
    torch.linalg.det = _patch_det
    torch.det = _patch_det  # Alias
    torch.linalg.cholesky = _patch_cholesky
    torch.linalg.eig = _patch_eig
    torch.linalg.eigh = _patch_eigh
    torch.linalg.lu_factor = _patch_lu_factor
    torch.linalg.lu_solve = _patch_lu_solve

    # Fix numpy compatibility
    torch.Tensor.__array__ = _tensor_array_compat
    torch.Tensor.numpy = lambda self: np.array(self.detach().tolist())

    # torch.from_numpy — WASM PyTorch build lacks NumPy interop
    _orig_from_numpy = torch.from_numpy

    def _patched_from_numpy(ndarray):
        try:
            return _orig_from_numpy(ndarray)
        except RuntimeError:
            return torch.tensor(ndarray.tolist())

    torch.from_numpy = _patched_from_numpy

    # torch.as_tensor — also uses numpy interop internally
    _orig_as_tensor = torch.as_tensor

    def _patched_as_tensor(data, dtype=None, device=None):
        try:
            return _orig_as_tensor(data, dtype=dtype, device=device)
        except (RuntimeError, TypeError):
            if hasattr(data, "tolist"):
                return torch.tensor(data.tolist(), dtype=dtype, device=device)
            return torch.tensor(data, dtype=dtype, device=device)

    torch.as_tensor = _patched_as_tensor

    # Tensor indexing — WASM PyTorch can't infer numpy dtypes for boolean/integer masks
    _orig_getitem = torch.Tensor.__getitem__

    def _patched_getitem(self, key):
        if isinstance(key, np.ndarray):
            key = torch.tensor(key.tolist())
        return _orig_getitem(self, key)

    torch.Tensor.__getitem__ = _patched_getitem

    _orig_setitem = torch.Tensor.__setitem__

    def _patched_setitem(self, key, value):
        if isinstance(key, np.ndarray):
            key = torch.tensor(key.tolist())
        return _orig_setitem(self, key, value)

    torch.Tensor.__setitem__ = _patched_setitem

    # torch.tensor — handle lists of 0-d tensors (WASM calls len() on elements, fails for 0-d)
    _orig_torch_tensor = torch.tensor

    def _patched_torch_tensor(data, *args, **kwargs):
        if isinstance(data, (list, tuple)):

            def _unwrap(item):
                if isinstance(item, torch.Tensor) and item.ndim == 0:
                    return item.item()
                return item

            data = type(data)(_unwrap(x) for x in data)
        return _orig_torch_tensor(data, *args, **kwargs)

    torch.tensor = _patched_torch_tensor

    print("✓ Torch linalg + numpy interop patches applied")


# ==============================================================================
# Torch testing patches
# ==============================================================================


class _LoggingTensorModeStub:
    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False


def _capture_logs_stub(*args, **kwargs):
    return _LoggingTensorModeStub()


def patch_torch_testing():
    """
    Stub torch.testing._internal modules that are missing in Pyodide.

    These modules are used by MACE imports but not needed for inference.
    """
    _internal = types.ModuleType("torch.testing._internal")
    _internal.__path__ = []
    _internal.__package__ = "torch.testing._internal"

    _common_utils = types.ModuleType("torch.testing._internal.common_utils")
    _common_utils.dtype_abbrs = {}

    _logging_tensor = types.ModuleType("torch.testing._internal.logging_tensor")
    _logging_tensor.LoggingTensorMode = _LoggingTensorModeStub
    _logging_tensor.capture_logs = _capture_logs_stub

    _internal.common_utils = _common_utils
    _internal.logging_tensor = _logging_tensor
    sys.modules["torch.testing._internal"] = _internal
    sys.modules["torch.testing._internal.common_utils"] = _common_utils
    sys.modules["torch.testing._internal.logging_tensor"] = _logging_tensor

    # Import torch.utils.checkpoint AFTER logging_tensor stubs are set up
    # (checkpoint.py imports logging_tensor at import time)
    import torch.utils.checkpoint  # noqa: F401

    print("✓ Torch testing patches applied")


# ==============================================================================
# Torch compiler patches
# ==============================================================================


def patch_torch_compiler():
    """
    Patch torch.compiler for Pyodide WASM.

    - Stubs torch.compiler.is_compiling (not available in WASM build)
    - Makes torch.compiler.disable a no-op decorator (avoids torch._dynamo
      import chain which requires C extensions missing in WASM)
    - Needed by e3nn >= 0.5 and fairchem-core which use @torch.compiler.disable
    """
    if not hasattr(torch, "compiler"):
        torch.compiler = types.ModuleType("torch.compiler")
    if not hasattr(torch.compiler, "is_compiling"):
        torch.compiler.is_compiling = lambda: False

    def _compiler_disable(fn=None, recursive=True):
        if fn is not None:
            return fn
        return lambda f: f

    torch.compiler.disable = _compiler_disable

    print("✓ Torch compiler patches applied")


# ==============================================================================
# Torch distributed patches (for FAIRChem)
# ==============================================================================


def patch_torch_distributed():
    """
    Stub torch.distributed modules that require C extensions missing in WASM.

    FAIRChem's mlip_unit.py imports from torch.distributed.checkpoint,
    torch.distributed.fsdp, and torch.distributed.device_mesh.
    All of these trigger torch._C._distributed_c10d which doesn't exist in WASM.
    """
    import enum

    import torch.distributed as _dist

    # Core distributed API stubs
    if not hasattr(_dist, "group"):

        class _Group:
            WORLD = None

        _dist.group = _Group
        _dist.WORLD = None
    if not hasattr(_dist, "ReduceOp"):

        class _ReduceOp:
            SUM = 0
            PRODUCT = 1
            MIN = 2
            MAX = 3
            BAND = 4
            BOR = 5
            BXOR = 6
            AVG = 7

        _dist.ReduceOp = _ReduceOp
    if not hasattr(_dist, "is_initialized"):
        _dist.is_initialized = lambda: False
    if not hasattr(_dist, "get_rank"):
        _dist.get_rank = lambda group=None: 0
    if not hasattr(_dist, "get_world_size"):
        _dist.get_world_size = lambda group=None: 1

    # torch.distributed.nn.functional
    _dist_nn = types.ModuleType("torch.distributed.nn")
    _dist_nn.__path__ = []
    _dist_nn.__package__ = "torch.distributed.nn"
    _dist_nn_func = types.ModuleType("torch.distributed.nn.functional")
    _dist_nn_func.__package__ = "torch.distributed.nn"
    _dist_nn_func.all_reduce = lambda tensor, *a, **k: tensor
    _dist_nn_func.reduce_scatter = lambda output, input_list, *a, **k: output
    _dist_nn_func.all_gather = lambda tensor_list, tensor, *a, **k: tensor_list
    _dist_nn.functional = _dist_nn_func
    sys.modules["torch.distributed.nn"] = _dist_nn
    sys.modules["torch.distributed.nn.functional"] = _dist_nn_func

    # C extension stubs
    sys.modules["torch._C._distributed_c10d"] = types.ModuleType("torch._C._distributed_c10d")
    _dc10d = types.ModuleType("torch.distributed.distributed_c10d")
    _dc10d.__package__ = "torch.distributed"
    sys.modules["torch.distributed.distributed_c10d"] = _dc10d

    # torch.distributed._shard
    _shard = types.ModuleType("torch.distributed._shard")
    _shard.__path__ = []
    _shard.__package__ = "torch.distributed._shard"
    sys.modules["torch.distributed._shard"] = _shard
    sys.modules["torch.distributed._shard.api"] = types.ModuleType("torch.distributed._shard.api")
    _shard_st = types.ModuleType("torch.distributed._shard.sharded_tensor")
    _shard_st.__path__ = []
    sys.modules["torch.distributed._shard.sharded_tensor"] = _shard_st
    _shard_st_meta = types.ModuleType("torch.distributed._shard.sharded_tensor.metadata")

    class _TensorProperties:
        pass

    _shard_st_meta.TensorProperties = _TensorProperties
    sys.modules["torch.distributed._shard.sharded_tensor.metadata"] = _shard_st_meta

    # torch.distributed.checkpoint
    _dcp = types.ModuleType("torch.distributed.checkpoint")
    _dcp.__path__ = []
    _dcp.__package__ = "torch.distributed.checkpoint"
    _dcp.save = lambda *a, **k: None
    _dcp.load = lambda *a, **k: None
    _dcp_meta = types.ModuleType("torch.distributed.checkpoint.metadata")
    _dcp_meta.TensorProperties = _TensorProperties

    class _BytesStorageMetadata:
        pass

    class _TensorStorageMetadata:
        pass

    class _Metadata:
        pass

    _dcp_meta.BytesStorageMetadata = _BytesStorageMetadata
    _dcp_meta.TensorStorageMetadata = _TensorStorageMetadata
    _dcp_meta.Metadata = _Metadata
    _dcp.metadata = _dcp_meta
    sys.modules["torch.distributed.checkpoint"] = _dcp
    sys.modules["torch.distributed.checkpoint.metadata"] = _dcp_meta

    for sub in [
        "state_dict",
        "stateful",
        "planner",
        "storage",
        "default_planner",
        "filesystem",
        "optimizer",
        "format_utils",
    ]:
        _sub_mod = types.ModuleType(f"torch.distributed.checkpoint.{sub}")
        _sub_mod.__package__ = "torch.distributed.checkpoint"
        sys.modules[f"torch.distributed.checkpoint.{sub}"] = _sub_mod
        setattr(_dcp, sub, _sub_mod)

    # format_utils stubs
    sys.modules["torch.distributed.checkpoint.format_utils"].dcp_to_torch_save = lambda *a, **k: None
    sys.modules["torch.distributed.checkpoint.format_utils"].torch_save_to_dcp = lambda *a, **k: None

    # state_dict stubs
    _sd_mod = sys.modules["torch.distributed.checkpoint.state_dict"]
    _sd_mod.get_model_state_dict = lambda model, *a, **k: model.state_dict() if hasattr(model, "state_dict") else {}
    _sd_mod.set_model_state_dict = (
        lambda model, sd, *a, **k: model.load_state_dict(sd) if hasattr(model, "load_state_dict") else None
    )
    _sd_mod.get_optimizer_state_dict = (
        lambda model, optim, *a, **k: optim.state_dict() if hasattr(optim, "state_dict") else {}
    )
    _sd_mod.get_state_dict = lambda model, *a, **k: model.state_dict() if hasattr(model, "state_dict") else {}
    _sd_mod.set_state_dict = lambda model, sd, *a, **k: None

    class _StateDictOptions:
        def __init__(self, **k):
            self.__dict__.update(k)

    _sd_mod.StateDictOptions = _StateDictOptions

    # stateful stub
    class _Stateful:
        pass

    sys.modules["torch.distributed.checkpoint.stateful"].Stateful = _Stateful

    # torch.distributed.fsdp
    _fsdp = types.ModuleType("torch.distributed.fsdp")
    _fsdp.__path__ = []
    _fsdp.__package__ = "torch.distributed.fsdp"
    sys.modules["torch.distributed.fsdp"] = _fsdp
    for fsdp_sub in ["fully_sharded_data_parallel", "api", "wrap", "sharded_grad_scaler"]:
        _fsub = types.ModuleType(f"torch.distributed.fsdp.{fsdp_sub}")
        _fsub.__package__ = "torch.distributed.fsdp"
        sys.modules[f"torch.distributed.fsdp.{fsdp_sub}"] = _fsub
        setattr(_fsdp, fsdp_sub, _fsub)

    # FSDP wrap policy stubs
    class _ModuleWrapPolicy:
        def __init__(self, module_classes=None):
            self.module_classes = module_classes or set()

    _fsdp_wrap = sys.modules["torch.distributed.fsdp.wrap"]
    _fsdp_wrap.ModuleWrapPolicy = _ModuleWrapPolicy
    _fsdp_wrap.lambda_auto_wrap_policy = lambda *a, **k: None
    _fsdp_wrap.transformer_auto_wrap_policy = lambda *a, **k: None

    # FSDP classes
    class _FullyShardedDataParallel(torch.nn.Module):
        def __init__(self, module, *a, **k):
            super().__init__()
            self.module = module

    class _ShardingStrategy(enum.Enum):
        FULL_SHARD = "FULL_SHARD"
        SHARD_GRAD_OP = "SHARD_GRAD_OP"
        NO_SHARD = "NO_SHARD"
        HYBRID_SHARD = "HYBRID_SHARD"

    class _MixedPrecision:
        def __init__(self, *a, **k):
            pass

    class _CPUOffload:
        def __init__(self, offload_params=False):
            self.offload_params = offload_params

    class _BackwardPrefetch(enum.Enum):
        BACKWARD_PRE = "BACKWARD_PRE"
        BACKWARD_POST = "BACKWARD_POST"

    class _StateDictTypeFSDP(enum.Enum):
        FULL_STATE_DICT = 0
        LOCAL_STATE_DICT = 1
        SHARDED_STATE_DICT = 2

    _fsdp.FullyShardedDataParallel = _FullyShardedDataParallel
    _fsdp.ShardingStrategy = _ShardingStrategy
    _fsdp.MixedPrecision = _MixedPrecision
    _fsdp.CPUOffload = _CPUOffload
    _fsdp.BackwardPrefetch = _BackwardPrefetch
    _fsdp.StateDictType = _StateDictTypeFSDP

    # FSDP StateDictConfig stubs
    class _StateDictConfig:
        def __init__(self, **k):
            self.__dict__.update(k)

    class _ShardedStateDictConfig(_StateDictConfig):
        pass

    class _FullStateDictConfig(_StateDictConfig):
        def __init__(self, offload_to_cpu=False, rank0_only=False, **k):
            super().__init__(**k)
            self.offload_to_cpu = offload_to_cpu
            self.rank0_only = rank0_only

    class _FullOptimStateDictConfig(_StateDictConfig):
        def __init__(self, offload_to_cpu=False, rank0_only=False, **k):
            super().__init__(**k)

    _fsdp.StateDictConfig = _StateDictConfig
    _fsdp.ShardedStateDictConfig = _ShardedStateDictConfig
    _fsdp.FullStateDictConfig = _FullStateDictConfig
    _fsdp.FullOptimStateDictConfig = _FullOptimStateDictConfig
    _fsdp.LocalStateDictConfig = type("LocalStateDictConfig", (_StateDictConfig,), {})
    _fsdp.OptimStateDictConfig = type("OptimStateDictConfig", (_StateDictConfig,), {})
    _fsdp.ShardedOptimStateDictConfig = type("ShardedOptimStateDictConfig", (_StateDictConfig,), {})

    # torch.distributed.device_mesh
    _dm = types.ModuleType("torch.distributed.device_mesh")
    _dm.__package__ = "torch.distributed"

    class _DeviceMesh:
        def __init__(self, *a, **k):
            pass

    _dm.DeviceMesh = _DeviceMesh
    _dm.init_device_mesh = lambda *a, **k: _DeviceMesh()
    sys.modules["torch.distributed.device_mesh"] = _dm

    # torch.distributed.tensor (DTensor)
    _dtensor = types.ModuleType("torch.distributed.tensor")
    _dtensor.__path__ = []
    _dtensor.__package__ = "torch.distributed.tensor"

    class _DTensor:
        pass

    _dtensor.DTensor = _DTensor
    sys.modules["torch.distributed.tensor"] = _dtensor

    # torch.distributed.algorithms
    _dalgo = types.ModuleType("torch.distributed.algorithms")
    _dalgo.__path__ = []
    _dalgo.__package__ = "torch.distributed.algorithms"
    sys.modules["torch.distributed.algorithms"] = _dalgo

    print("✓ Torch distributed patches applied")


# ==============================================================================
# FAIRChem heavy dependency stubs
# ==============================================================================


def _make_stub_module(name, attrs=None, submodules=None):
    """Create a stub module with optional attributes and submodules."""
    from importlib.machinery import ModuleSpec

    mod = types.ModuleType(name)
    mod.__path__ = []
    mod.__package__ = name
    mod.__version__ = "0.0.0"
    mod.__spec__ = ModuleSpec(name, None, is_package=True)
    if attrs:
        for k, v in attrs.items():
            setattr(mod, k, v)
    sys.modules[name] = mod
    if submodules:
        for sub_name in submodules:
            full_name = f"{name}.{sub_name}"
            sub_mod = types.ModuleType(full_name)
            sub_mod.__path__ = []
            sub_mod.__package__ = name
            sub_mod.__spec__ = ModuleSpec(full_name, None, is_package=True)
            setattr(mod, sub_name, sub_mod)
            sys.modules[full_name] = sub_mod
    return mod


def patch_fairchem_deps():
    """
    Stub heavy dependencies that fairchem-core imports but doesn't need for inference.

    This stubs: numba, ray (+ serve), wandb, torchtnt, hydra, omegaconf,
    submitit, clusterscope, tqdm, huggingface_hub, websockets.
    """
    # --- numba ---
    numba_mod = _make_stub_module("numba", submodules=["core", "core.types", "typed"])
    numba_mod.njit = lambda *a, **k: (lambda f: f) if not a or callable(a[0]) else lambda f: f
    numba_mod.jit = numba_mod.njit
    numba_mod.prange = range
    for t in ("int32", "int64", "float32", "float64", "boolean"):
        setattr(numba_mod, t, t)

    class _TypedList(list):
        pass

    sys.modules["numba.typed"].List = _TypedList

    # --- ray ---
    ray_mod = _make_stub_module(
        "ray",
        submodules=[
            "serve",
            "runtime_env",
            "train",
            "data",
            "util",
            "util.scheduling_strategies",
            "util.queue",
        ],
    )

    def _ray_remote(*args, **kwargs):
        if args and callable(args[0]):
            args[0]._remote = lambda *a, **kw: None
            return args[0]

        def wrapper(fn_or_cls):
            fn_or_cls._remote = lambda *a, **kw: None
            return fn_or_cls

        return wrapper

    ray_mod.remote = _ray_remote
    ray_mod.init = lambda *a, **k: None
    ray_mod.get = lambda *a, **k: None
    ray_mod.put = lambda *a, **k: None
    ray_mod.wait = lambda *a, **k: ([], [])
    ray_mod.is_initialized = lambda: False

    class _ObjectRef:
        pass

    ray_mod.ObjectRef = _ObjectRef

    class _PlacementGroupSchedulingStrategy:
        def __init__(self, *a, **k):
            pass

    sys.modules["ray.util.scheduling_strategies"].PlacementGroupSchedulingStrategy = (
        _PlacementGroupSchedulingStrategy
    )

    # ray.serve stubs
    _serve = sys.modules["ray.serve"]

    def _serve_deployment(*args, **kwargs):
        if args and callable(args[0]):
            return args[0]
        return lambda cls_or_fn: cls_or_fn

    _serve.deployment = _serve_deployment
    _serve.ingress = lambda *a, **k: (lambda cls: cls)
    _serve.run = lambda *a, **k: None
    _serve.batch = lambda *args, **kwargs: (lambda fn: fn) if not args or not callable(args[0]) else args[0]

    _serve_schema = types.ModuleType("ray.serve.schema")
    _serve_schema.__package__ = "ray.serve"

    class _LoggingConfig:
        def __init__(self, **k):
            self.__dict__.update(k)

    _serve_schema.LoggingConfig = _LoggingConfig
    _serve.schema = _serve_schema
    sys.modules["ray.serve.schema"] = _serve_schema

    # --- wandb ---
    wandb_mod = _make_stub_module("wandb")
    wandb_mod.init = lambda *a, **k: None
    wandb_mod.log = lambda *a, **k: None
    wandb_mod.finish = lambda *a, **k: None

    # --- torchtnt ---
    _make_stub_module(
        "torchtnt",
        submodules=[
            "framework",
            "framework.state",
            "framework.unit",
            "framework.callback",
            "framework.auto_unit",
            "framework.fit",
            "framework.train",
            "framework.evaluate",
            "framework.predict",
            "utils",
            "utils.loggers",
            "utils.timer",
            "utils.distributed",
            "utils.prepare_module",
        ],
    )

    class _PredictUnit:
        def __init__(self, *a, **k):
            pass

        def __class_getitem__(cls, item):
            return cls

    class _TrainUnit:
        def __init__(self, *a, **k):
            pass

        def __class_getitem__(cls, item):
            return cls

    class _EvalUnit:
        def __init__(self, *a, **k):
            pass

        def __class_getitem__(cls, item):
            return cls

    class _State:
        def __init__(self, *a, **k):
            pass

        def __class_getitem__(cls, item):
            return cls

    class _Callback:
        pass

    tnt_framework = sys.modules["torchtnt.framework"]
    tnt_unit = sys.modules["torchtnt.framework.unit"]
    tnt_state = sys.modules["torchtnt.framework.state"]
    tnt_cb = sys.modules["torchtnt.framework.callback"]
    for m in [tnt_framework, tnt_unit]:
        m.PredictUnit = _PredictUnit
        m.TrainUnit = _TrainUnit
        m.EvalUnit = _EvalUnit
    tnt_state.State = _State
    tnt_cb.Callback = _Callback
    tnt_framework.State = _State
    tnt_framework.Callback = _Callback

    # torchtnt entry point functions
    sys.modules["torchtnt.framework.fit"].fit = lambda *a, **k: None
    sys.modules["torchtnt.framework.train"].train = lambda *a, **k: None
    sys.modules["torchtnt.framework.evaluate"].evaluate = lambda *a, **k: None
    sys.modules["torchtnt.framework.predict"].predict = lambda *a, **k: None

    # torchtnt.utils stubs
    tnt_dist = sys.modules["torchtnt.utils.distributed"]
    tnt_dist.get_file_init_method = lambda *a, **k: ""
    tnt_dist.get_tcp_init_method = lambda *a, **k: ""
    tnt_dist.spawn_multi_process = lambda *a, **k: None

    tnt_prep = sys.modules["torchtnt.utils.prepare_module"]
    tnt_prep.prepare_module = lambda module, *a, **k: module
    tnt_prep.FSDPStrategy = type("FSDPStrategy", (), {"__init__": lambda self, **k: None})
    tnt_prep.DDPStrategy = type("DDPStrategy", (), {"__init__": lambda self, **k: None})
    tnt_prep.NOOPStrategy = type("NOOPStrategy", (), {"__init__": lambda self, **k: None})

    # --- hydra / omegaconf ---
    omegaconf_mod = _make_stub_module("omegaconf")

    class _DictConfig(dict):
        pass

    class _ListConfig(list):
        pass

    omegaconf_mod.DictConfig = _DictConfig
    omegaconf_mod.ListConfig = _ListConfig
    omegaconf_mod.OmegaConf = type(
        "OmegaConf",
        (),
        {
            "to_container": staticmethod(lambda cfg, **k: dict(cfg) if isinstance(cfg, dict) else cfg),
            "create": staticmethod(lambda d: _DictConfig(d) if isinstance(d, dict) else d),
        },
    )
    _make_stub_module("hydra", submodules=["core", "core.global_hydra", "utils"])

    def _hydra_instantiate(config, *args, _recursive_=True, **kwargs):
        import importlib
        if isinstance(config, dict) and "_target_" in config:
            target = config["_target_"]
            mod_path, cls_name = target.rsplit(".", 1)
            mod = importlib.import_module(mod_path)
            cls = getattr(mod, cls_name)
            pos_args = list(args) + list(config.get("_args_", []))
            cfg = {k: v for k, v in config.items() if not k.startswith("_")}
            if _recursive_:
                for k, v in cfg.items():
                    if isinstance(v, dict) and "_target_" in v:
                        cfg[k] = _hydra_instantiate(v)
                    elif isinstance(v, list):
                        cfg[k] = [_hydra_instantiate(i) if isinstance(i, dict) and "_target_" in i else i for i in v]
                pos_args = [_hydra_instantiate(a) if isinstance(a, dict) and "_target_" in a else a for a in pos_args]
            return cls(*pos_args, **{**cfg, **kwargs})
        return config

    sys.modules["hydra.utils"].instantiate = _hydra_instantiate
    sys.modules["hydra"].utils = sys.modules["hydra.utils"]

    # --- submitit / clusterscope ---
    _make_stub_module("submitit")
    _make_stub_module("clusterscope")

    # --- websockets ---
    _make_stub_module("websockets")

    # --- tqdm ---
    tqdm_mod = _make_stub_module("tqdm", submodules=["auto", "std"])

    def _tqdm_passthrough(iterable=None, *a, **k):
        return iterable if iterable is not None else iter([])

    tqdm_mod.tqdm = _tqdm_passthrough
    sys.modules["tqdm.auto"].tqdm = _tqdm_passthrough
    sys.modules["tqdm.std"].tqdm = _tqdm_passthrough

    # --- huggingface_hub ---
    hf_mod = _make_stub_module("huggingface_hub", submodules=["utils"])
    hf_mod.hf_hub_download = lambda *a, **k: ""
    hf_mod.snapshot_download = lambda *a, **k: ""

    # --- ase_db_backends ---
    _make_stub_module("ase_db_backends")
    # --- INT8 quantized model support ---
    _orig_torch_load = torch.load

    def _int8_aware_torch_load(f, *args, **kwargs):
        result = _orig_torch_load(f, *args, **kwargs)
        if isinstance(result, dict) and "quantized_ema_state_dict" in result:
            import gc as _gc
            from fairchem.core.units.mlip_unit.api.inference import MLIPInferenceCheckpoint

            print("  Dequantizing INT8 → FP16 (streaming)...")
            quantized_ema = result.pop("quantized_ema_state_dict")
            scales = result.pop("quantization_scales")
            ema_state_dict = {}
            names = list(quantized_ema.keys())
            for name in names:
                tensor = quantized_ema.pop(name)
                if name in scales:
                    scale = scales.pop(name)
                    ema_state_dict[name] = (tensor.float() * scale.float()).half()
                    del scale
                else:
                    ema_state_dict[name] = tensor
                del tensor
            del quantized_ema, scales
            _gc.collect()
            checkpoint = MLIPInferenceCheckpoint(
                model_config=result["model_config"],
                model_state_dict=result.get("model_state_dict", {}),
                ema_state_dict=ema_state_dict,
                tasks_config=result["tasks_config"],
            )
            del result
            _gc.collect()
            print("  ✓ Dequantization complete")
            return checkpoint
        return result

    torch.load = _int8_aware_torch_load

    # --- Model registry fallback ---
    try:
        from fairchem.core.common.registry import registry as _registry

        _orig_get_model = _registry.get_model_class

        def _fallback_get_model_class(name):
            try:
                return _orig_get_model(name)
            except RuntimeError:
                import importlib as _imp
                module_path, class_name = name.rsplit(".", 1)
                mod = _imp.import_module(module_path)
                return getattr(mod, class_name)

        _registry.get_model_class = _fallback_get_model_class
    except Exception:
        pass

    print("✓ FAIRChem dependency stubs applied")


# ==============================================================================
# Matscipy patches
# ==============================================================================


def _matscipy_neighbour_list_compat(quantities, atoms=None, cutoff=None, positions=None, cell=None, pbc=None, **_):
    """Translate matscipy-style keyword-arg call into an ASE neighbor_list call."""
    from ase import Atoms as _Atoms
    from ase.neighborlist import neighbor_list as _ase_neighbor_list

    if atoms is None:
        atoms = _Atoms(positions=positions, cell=cell, pbc=pbc if pbc is not None else [False, False, False])
    return _ase_neighbor_list(quantities, atoms, cutoff)


def patch_matscipy():
    """
    Stub matscipy package with ASE neighbor_list compatibility.

    Matscipy is a C-extension package that cannot be compiled for WASM.
    MACE uses it for neighbor calculations, which we redirect to ASE.
    """
    _matscipy = types.ModuleType("matscipy")
    _matscipy.__path__ = []
    _matscipy.__package__ = "matscipy"
    _matscipy_neighbours = types.ModuleType("matscipy.neighbours")
    _matscipy_neighbours.neighbour_list = _matscipy_neighbour_list_compat
    _matscipy.neighbours = _matscipy_neighbours
    sys.modules["matscipy"] = _matscipy
    sys.modules["matscipy.neighbours"] = _matscipy_neighbours

    print("✓ Matscipy patches applied")


# ==============================================================================
# LMDB and HDF5 patches
# ==============================================================================


def patch_mace_training():
    """
    Stub lmdb and h5py packages.

    These are C-extension packages used by MACE's training/dataset code
    but not needed for inference. Stubs allow imports to succeed.
    """
    for _pkg in ("lmdb", "h5py"):
        if _pkg not in sys.modules:
            sys.modules[_pkg] = types.ModuleType(_pkg)

    print("✓ LMDB and HDF5 stubs applied")


# ==============================================================================
# MACE tools patches
# ==============================================================================


def patch_mace_tools():
    """
    Fix MACE's torch_geometric import order issues in Pyodide.

    In Pyodide, torch_geometric.data may not be set during circular imports.
    Pre-importing ensures the attribute is available when MACE needs it.
    """
    try:
        import importlib as _importlib

        _tg = _importlib.import_module("mace.tools.torch_geometric")
        _tg_data = _importlib.import_module("mace.tools.torch_geometric.data")
        _tg.data = _tg_data
        print("✓ MACE tools patches applied")
    except Exception as e:
        print(f"⚠ MACE tools patches skipped: {e}")


# ==============================================================================
# Convenience function to apply all patches
# ==============================================================================


def apply_all_patches(include_fairchem=False, include_mattersim=False, include_sevennet=False):
    """
    Apply all torch and model patches for Pyodide in one call.

    Args:
        include_fairchem: If True, also apply FAIRChem-specific patches
            (torch.distributed, heavy dependency stubs). Set this when
            using fairchem-core / UMA models.
        include_mattersim: If True, also apply MatterSim-specific patches
            (loguru, azure, e3nn JIT stubs). Set this when
            using MatterSim / M3GNet models.
        include_sevennet: If True, also apply SevenNet-specific patches
            (pandas, tqdm, packaging stubs). Set this when
            using SevenNet / 7net models.
    """
    patch_torch_linalg()
    patch_torch_compiler()
    patch_torch_testing()
    patch_matscipy()
    patch_mace_training()
    patch_mace_tools()

    if include_fairchem:
        patch_torch_distributed()
        patch_fairchem_deps()

    if include_mattersim:
        patch_torch_distributed()
        patch_mattersim_deps()

    if include_sevennet:
        patch_torch_distributed()
        patch_sevennet_deps()

    print("\n✅ All Pyodide patches applied successfully!")


# ==============================================================================
# MatterSim patches
# ==============================================================================


def patch_mattersim_deps():
    """
    Stub heavy dependencies required by MatterSim but not needed for inference.

    Stubs: loguru, azure.*, atomate2, seekpath, phonopy, phono3py, mp_api,
    sklearn, and patches e3nn to disable JIT/torch.compile.
    """
    # --- loguru ---
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

    # --- azure (not needed for inference) ---
    _make_stub_module("azure", submodules=[
        "identity", "storage", "storage.blob",
    ])

    # --- heavy optional deps (not needed for inference) ---
    for pkg in [
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
        _make_stub_module(pkg)

    # --- scikit-learn (stub) ---
    sk_mod = _make_stub_module("sklearn", submodules=[
        "base", "utils", "utils.validation",
        "preprocessing", "model_selection",
        "gaussian_process", "gaussian_process.kernels",
    ])

    # Stub GaussianProcessRegressor
    class _GPR:
        def __init__(self, *a, **k):
            pass
        def fit(self, *a, **k):
            return self
        def predict(self, X, return_std=False):
            import numpy as _np
            mean = _np.zeros(X.shape[0])
            if return_std:
                return mean, _np.ones(X.shape[0])
            return mean
        def log_marginal_likelihood(self):
            return 0.0
    sys.modules["sklearn.gaussian_process"].GaussianProcessRegressor = _GPR

    # Stub kernels
    class _Kernel:
        pass
    class _DotProduct(_Kernel):
        def __init__(self, *a, **k):
            pass
    class _Hyperparameter:
        def __init__(self, *a, **k):
            pass
    _sk_kernels = sys.modules["sklearn.gaussian_process.kernels"]
    _sk_kernels.Kernel = _Kernel
    _sk_kernels.DotProduct = _DotProduct
    _sk_kernels.Hyperparameter = _Hyperparameter

    # --- requests (for pyodide-http compat) ---
    try:
        import pyodide_http  # noqa: F401
        pyodide_http.patch_all()
    except ImportError:
        pass

    # --- patch e3nn to skip JIT compilation ---
    try:
        import e3nn
        e3nn._SO3_INITIALIZED = True  # skip init
    except Exception:
        pass

    # --- patch torch.jit for e3nn ---
    import torch
    if not hasattr(torch.jit, '_original_script'):
        _orig_script = torch.jit.script

        def _noop_script(obj=None, *a, **k):
            if obj is not None:
                return obj
            return lambda fn: fn

        torch.jit.script = _noop_script

    # --- torch_ema (training only) ---
    _te = _make_stub_module("torch_ema")

    class _EMA:
        def __init__(self, *a, **k):
            pass
    _te.ExponentialMovingAverage = _EMA

    # --- torchmetrics (training only) ---
    _tm = _make_stub_module("torchmetrics")

    class _MeanMetric:
        def __init__(self, *a, **k):
            pass
    _tm.MeanMetric = _MeanMetric

    # --- torch_geometric (data loading only) ---
    _tg = _make_stub_module("torch_geometric", submodules=[
        "data", "loader", "utils",
    ])

    class _Data:
        def __init__(self, **kwargs):
            for k, v in kwargs.items():
                setattr(self, k, v)

        def to(self, device):
            import torch as _t
            for k, v in self.__dict__.items():
                if isinstance(v, _t.Tensor):
                    setattr(self, k, v.to(device))
            return self

    sys.modules["torch_geometric.data"].Data = _Data

    class _DataLoader:
        """Minimal DataLoader that yields items with Batch-style attributes."""
        def __init__(self, dataset, batch_size=1, shuffle=False, **kwargs):
            self._dataset = list(dataset)

        def __iter__(self):
            import torch as _t
            for item in self._dataset:
                # Convert scalar int/float fields to 1-element tensors
                # (torch_geometric Batch collation does this automatically)
                for attr in list(vars(item).keys()):
                    val = getattr(item, attr)
                    if isinstance(val, (int, float)):
                        setattr(item, attr, _t.tensor([val]))
                # Add torch_geometric Batch-style attributes for single graph
                if not hasattr(item, 'num_graphs'):
                    item.num_graphs = 1
                if not hasattr(item, 'batch'):
                    n_atoms = item.num_atoms if hasattr(item, 'num_atoms') else _t.tensor([0])
                    if isinstance(n_atoms, _t.Tensor):
                        n_atoms = int(n_atoms.item())
                    item.batch = _t.zeros(n_atoms, dtype=_t.long)
                yield item

        def __len__(self):
            return len(self._dataset)

    sys.modules["torch_geometric.loader"].DataLoader = _DataLoader

    # --- torch_runstats (scatter ops) ---
    _trs = _make_stub_module("torch_runstats", submodules=["scatter"])
    import torch as _torch

    def _scatter(src, index, dim_size=None, dim=0, reduce="sum"):
        if dim_size is None:
            dim_size = int(index.max()) + 1
        out = _torch.zeros(dim_size, *src.shape[1:], dtype=src.dtype, device=src.device)
        if src.dim() == 1:
            idx = index
        else:
            idx = index.unsqueeze(-1).expand_as(src)
        if reduce == "sum" or reduce == "add":
            out.scatter_add_(0, idx, src)
        elif reduce == "mean":
            out.scatter_add_(0, idx, src)
            count = _torch.zeros(dim_size, dtype=src.dtype, device=src.device)
            count.scatter_add_(0, index, _torch.ones(index.shape[0], dtype=src.dtype, device=src.device))
            count = count.clamp(min=1)
            if src.dim() > 1:
                count = count.unsqueeze(-1)
            out = out / count
        return out

    sys.modules["torch_runstats.scatter"].scatter = _scatter

    def _scatter_mean(src, index, dim_size=None, dim=0):
        return _scatter(src, index, dim_size=dim_size, dim=dim, reduce="mean")

    sys.modules["torch_runstats.scatter"].scatter_mean = _scatter_mean

    print("✓ MatterSim dependency stubs applied")


# ==============================================================================
# SevenNet patches
# ==============================================================================


def patch_sevennet_deps():
    """
    Stub heavy dependencies required by SevenNet but not needed for inference.

    Stubs: pandas, braceexpand, tqdm, packaging, and patches
    torch_geometric.data.Data for SevenNet's AtomGraphData.
    """
    import sys

    # --- packaging (used in __init__.py version check) ---
    # Only stub if real packaging is not available
    if "packaging.version" not in sys.modules:
        try:
            import packaging.version  # noqa: F401
        except (ImportError, ModuleNotFoundError):
            pkg = _make_stub_module("packaging", submodules=["version"])

            class _Version:
                def __init__(self, v):
                    self._v = str(v)
                    parts = self._v.split(".")
                    self.major = int(parts[0]) if parts else 0
                    self.minor = int(parts[1]) if len(parts) > 1 else 0
                    self.micro = int(parts[2]) if len(parts) > 2 else 0

                def __lt__(self, other):
                    return (self.major, self.minor, self.micro) < (other.major, other.minor, other.micro)

                def __ge__(self, other):
                    return not self.__lt__(other)

                def __repr__(self):
                    return f"Version('{self._v}')"

            sys.modules["packaging.version"].Version = _Version
            sys.modules["packaging.version"].parse = lambda v: _Version(v)

    # --- pandas (used in checkpoint.py, not needed for inference) ---
    try:
        import pandas  # noqa: F401
    except (ImportError, ModuleNotFoundError):
        pd = _make_stub_module("pandas", submodules=["core", "core.frame"])

        class _DataFrame:
            def __init__(self, *a, **k):
                pass
        pd.DataFrame = _DataFrame
        sys.modules["pandas.core.frame"].DataFrame = _DataFrame

    # --- braceexpand (used in train/dataload.py) ---
    try:
        import braceexpand  # noqa: F401
    except (ImportError, ModuleNotFoundError):
        be = _make_stub_module("braceexpand")
        be.braceexpand = lambda s: [s]

    # --- tqdm (used in util.py and train) ---
    try:
        import tqdm  # noqa: F401
    except (ImportError, ModuleNotFoundError):
        tqdm_mod = _make_stub_module("tqdm", submodules=["auto"])

        class _tqdm:
            def __init__(self, iterable=None, *a, **k):
                self._it = iterable

            def __iter__(self):
                return iter(self._it) if self._it else iter([])

            def __enter__(self):
                return self

            def __exit__(self, *a):
                pass

            def update(self, *a):
                pass

            def close(self):
                pass

        tqdm_mod.tqdm = _tqdm
        sys.modules["tqdm.auto"].tqdm = _tqdm

    # --- torch_geometric.data.Data needs proper item access for SevenNet ---
    import torch
    tg_data = sys.modules.get("torch_geometric.data")
    if tg_data is None:
        _make_stub_module("torch_geometric", submodules=["data", "loader", "utils"])
        tg_data = sys.modules["torch_geometric.data"]

    class _SevenNetData:
        """torch_geometric-compatible Data with dict-like access for SevenNet."""
        def __init__(self, x=None, edge_index=None, edge_attr=None, y=None, pos=None, **kwargs):
            self._store = {}
            if x is not None:
                self._store['x'] = x
            if edge_index is not None:
                self._store['edge_index'] = edge_index
            if edge_attr is not None:
                self._store['edge_attr'] = edge_attr
            if y is not None:
                self._store['y'] = y
            if pos is not None:
                self._store['pos'] = pos
            for k, v in kwargs.items():
                self._store[k] = v

        def __setitem__(self, key, value):
            self._store[key] = value

        def __getitem__(self, key):
            return self._store[key]

        def __contains__(self, key):
            return key in self._store

        def __delitem__(self, key):
            del self._store[key]

        def __setattr__(self, name, value):
            if name == '_store':
                super().__setattr__(name, value)
            else:
                self._store[name] = value

        def __getattr__(self, name):
            if name == '_store':
                raise AttributeError(name)
            try:
                return self._store[name]
            except KeyError:
                raise AttributeError(name)

        def to(self, device):
            for k, v in self._store.items():
                if isinstance(v, torch.Tensor):
                    self._store[k] = v.to(device)
            return self

        def to_dict(self):
            return dict(self._store)

        def keys(self):
            return self._store.keys()

        @classmethod
        def from_numpy_dict(cls, d):
            """Convert numpy arrays to tensors."""
            import numpy as np
            obj = cls()
            for k, v in d.items():
                if isinstance(v, np.ndarray):
                    if v.dtype in (np.float32, np.float64):
                        obj._store[k] = torch.from_numpy(v).float()
                    elif v.dtype in (np.int32, np.int64):
                        obj._store[k] = torch.from_numpy(v).long()
                    else:
                        obj._store[k] = torch.from_numpy(v)
                elif isinstance(v, (int, float)):
                    obj._store[k] = v
                else:
                    obj._store[k] = v
            return obj

    tg_data.Data = _SevenNetData

    # --- patch e3nn to skip JIT compilation ---
    try:
        import e3nn
        e3nn._SO3_INITIALIZED = True
        # Set eager mode to avoid torch.jit.script (not supported in Pyodide)
        if hasattr(e3nn, '_OPT_DEFAULTS'):
            e3nn._OPT_DEFAULTS["jit_mode"] = "eager"
    except Exception:
        pass

    # --- patch torch.jit for e3nn ---
    if not hasattr(torch.jit, '_original_script'):
        _orig_script = torch.jit.script

        def _noop_script(obj=None, *a, **k):
            if obj is not None:
                return obj
            return lambda fn: fn

        torch.jit.script = _noop_script

    # --- patch e3nn compile_mode decorator ---
    try:
        import e3nn.util.jit
        e3nn.util.jit.compile_mode = lambda mode: lambda cls: cls
    except Exception:
        pass

    print("✓ SevenNet dependency stubs applied")
