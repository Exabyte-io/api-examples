import enum
import gc
import sys
import types
from importlib import import_module

import torch

from ...primitive.environment import is_pyodide_environment
from .torch import _ensure_omegaconf_stub, _make_stub_module, patch_torch_distributed


def apply_patches():
    patch_torch_distributed()
    patch_fairchem_deps()


def patch_fairchem_deps():
    """
    Stub heavy dependencies that fairchem-core imports but doesn't need for inference.

    This stubs: numba, ray (+ serve), wandb, torchtnt, hydra, omegaconf,
    submitit, psutil, clusterscope, tqdm, huggingface_hub, websockets.
    """
    numba_mod = _make_stub_module("numba", submodules=["core", "core.types", "typed"])
    numba_mod.njit = lambda *a, **k: (lambda f: f) if not a or callable(a[0]) else lambda f: f
    numba_mod.jit = numba_mod.njit
    numba_mod.prange = range
    for t in ("int32", "int64", "float32", "float64", "boolean"):
        setattr(numba_mod, t, t)

    class _TypedList(list):
        pass

    sys.modules["numba.typed"].List = _TypedList
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

    sys.modules["ray.util.scheduling_strategies"].PlacementGroupSchedulingStrategy = _PlacementGroupSchedulingStrategy
    _patch_ray_serve()
    _patch_wandb()
    _patch_torchtnt()
    _patch_hydra()
    _patch_submitit()
    _patch_psutil()
    _make_stub_module("clusterscope")
    _make_stub_module("websockets")
    _patch_tqdm()
    _patch_huggingface_hub()
    _make_stub_module("ase_db_backends")
    _patch_int8_torch_load()
    _patch_model_registry_fallback()

    print("✓ FAIRChem dependency stubs applied")


def _patch_ray_serve():
    serve_mod = sys.modules["ray.serve"]

    def _serve_deployment(*args, **kwargs):
        if args and callable(args[0]):
            return args[0]
        return lambda cls_or_fn: cls_or_fn

    serve_mod.deployment = _serve_deployment
    serve_mod.ingress = lambda *a, **k: (lambda cls: cls)
    serve_mod.run = lambda *a, **k: None
    serve_mod.batch = lambda *args, **kwargs: (lambda fn: fn) if not args or not callable(args[0]) else args[0]
    serve_mod.multiplexed = lambda *args, **kwargs: (lambda fn: fn) if not args or not callable(args[0]) else args[0]

    serve_schema = types.ModuleType("ray.serve.schema")
    serve_schema.__package__ = "ray.serve"

    class _LoggingConfig:
        def __init__(self, **k):
            self.__dict__.update(k)

    class _ApplicationStatus(enum.Enum):
        NOT_STARTED = "NOT_STARTED"
        DEPLOYING = "DEPLOYING"
        RUNNING = "RUNNING"
        DEPLOY_FAILED = "DEPLOY_FAILED"
        DELETING = "DELETING"

    serve_schema.LoggingConfig = _LoggingConfig
    serve_schema.ApplicationStatus = _ApplicationStatus
    serve_mod.schema = serve_schema
    sys.modules["ray.serve.schema"] = serve_schema


def _patch_submitit():
    submitit_mod = _make_stub_module("submitit", submodules=["helpers", "core", "core.utils"])

    class _SubmititPlaceholder:
        pass

    submitit_mod.Job = _SubmititPlaceholder
    sys.modules["submitit.helpers"].Checkpointable = _SubmititPlaceholder
    sys.modules["submitit.helpers"].DelayedSubmission = _SubmititPlaceholder
    sys.modules["submitit.core.utils"].JobPaths = _SubmititPlaceholder


def _patch_psutil():
    psutil_mod = types.ModuleType("psutil")

    class _UnavailableProcess:
        def __init__(self, *args, **kwargs):
            raise RuntimeError("psutil process management is unavailable in Pyodide")

    psutil_mod.Process = _UnavailableProcess
    psutil_mod.wait_procs = lambda *args, **kwargs: ([], [])
    sys.modules["psutil"] = psutil_mod


def _patch_wandb():
    wandb_mod = _make_stub_module("wandb")
    wandb_mod.init = lambda *a, **k: None
    wandb_mod.log = lambda *a, **k: None
    wandb_mod.finish = lambda *a, **k: None


def _patch_torchtnt():
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
    for module in [tnt_framework, tnt_unit]:
        module.PredictUnit = _PredictUnit
        module.TrainUnit = _TrainUnit
        module.EvalUnit = _EvalUnit
    tnt_state.State = _State
    tnt_cb.Callback = _Callback
    tnt_framework.State = _State
    tnt_framework.Callback = _Callback
    sys.modules["torchtnt.framework.fit"].fit = lambda *a, **k: None
    sys.modules["torchtnt.framework.train"].train = lambda *a, **k: None
    sys.modules["torchtnt.framework.evaluate"].evaluate = lambda *a, **k: None
    sys.modules["torchtnt.framework.predict"].predict = lambda *a, **k: None

    tnt_dist = sys.modules["torchtnt.utils.distributed"]
    tnt_dist.get_file_init_method = lambda *a, **k: ""
    tnt_dist.get_tcp_init_method = lambda *a, **k: ""
    tnt_dist.spawn_multi_process = lambda *a, **k: None

    tnt_prep = sys.modules["torchtnt.utils.prepare_module"]
    tnt_prep.prepare_module = lambda module, *a, **k: module
    tnt_prep.FSDPStrategy = type("FSDPStrategy", (), {"__init__": lambda self, **k: None})
    tnt_prep.DDPStrategy = type("DDPStrategy", (), {"__init__": lambda self, **k: None})
    tnt_prep.NOOPStrategy = type("NOOPStrategy", (), {"__init__": lambda self, **k: None})


def _patch_hydra():
    _ensure_omegaconf_stub()
    _make_stub_module("hydra", submodules=["core", "core.global_hydra", "utils"])

    def _hydra_instantiate(config, *args, _recursive_=True, **kwargs):
        if isinstance(config, dict) and "_target_" in config:
            target = config["_target_"]
            mod_path, cls_name = target.rsplit(".", 1)
            mod = import_module(mod_path)
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


def _patch_tqdm():
    tqdm_mod = _make_stub_module("tqdm", submodules=["auto", "std"])

    def _tqdm_passthrough(iterable=None, *a, **k):
        return iterable if iterable is not None else iter([])

    tqdm_mod.tqdm = _tqdm_passthrough
    sys.modules["tqdm.auto"].tqdm = _tqdm_passthrough
    sys.modules["tqdm.std"].tqdm = _tqdm_passthrough


def _patch_huggingface_hub():
    hf_mod = _make_stub_module("huggingface_hub", submodules=["utils"])
    hf_mod.hf_hub_download = lambda *a, **k: ""
    hf_mod.snapshot_download = lambda *a, **k: ""


def _patch_int8_torch_load():
    original_torch_load = torch.load

    def _int8_aware_torch_load(f, *args, **kwargs):
        result = original_torch_load(f, *args, **kwargs)
        if isinstance(result, dict) and "quantized_ema_state_dict" in result:
            return _dequantize_int8_checkpoint(result)
        return result

    torch.load = _int8_aware_torch_load


def _dequantize_int8_checkpoint(result):
    from fairchem.core.units.mlip_unit.api.inference import MLIPInferenceCheckpoint

    print("  Dequantizing INT8 -> FP16 (streaming)...")
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
    gc.collect()
    checkpoint = MLIPInferenceCheckpoint(
        model_config=result["model_config"],
        model_state_dict=result.get("model_state_dict", {}),
        ema_state_dict=ema_state_dict,
        tasks_config=result["tasks_config"],
    )
    del result
    gc.collect()
    print("  ✓ Dequantization complete")
    return checkpoint


def _patch_model_registry_fallback():
    try:
        from fairchem.core.common.registry import registry as fairchem_registry

        original_get_model = fairchem_registry.get_model_class

        def _fallback_get_model_class(name):
            try:
                return original_get_model(name)
            except RuntimeError:
                module_path, class_name = name.rsplit(".", 1)
                mod = import_module(module_path)
                return getattr(mod, class_name)

        fairchem_registry.get_model_class = _fallback_get_model_class
    except Exception:
        pass


MODEL_PATHS_MAP = {
    "f16": "/drive/packages/models/uma-s-1p1-f16.pt",
    "int8": "/drive/packages/models/uma-s-1p1-int8.pt",
}


def get_uma_model_pyodide(model: str, task_name="omat", device="cpu", **kwargs):
    if model not in MODEL_PATHS_MAP:
        raise ValueError(f"Invalid model name: {model}. Valid options are: {list(MODEL_PATHS_MAP.keys())}")
    fairchem_core = import_module("fairchem.core")
    mlip_unit = import_module("fairchem.core.units.mlip_unit")
    predictor = mlip_unit.load_predict_unit(MODEL_PATHS_MAP[model], device=device, **kwargs)
    return fairchem_core.FAIRChemCalculator(predictor, task_name=task_name)


def create_uma_calculator(model="f16", task_name="omat", device="cpu", model_path=None, checkpoint=None, **kwargs):
    if is_pyodide_environment():
        return get_uma_model_pyodide(model=model, task_name=task_name, device=device, **kwargs)

    resolved_model_path = model_path or checkpoint

    fairchem_core = import_module("fairchem.core")
    mlip_unit = import_module("fairchem.core.units.mlip_unit")
    predictor = mlip_unit.load_predict_unit(str(resolved_model_path), device=device, **kwargs)
    return fairchem_core.FAIRChemCalculator(predictor, task_name=task_name)
