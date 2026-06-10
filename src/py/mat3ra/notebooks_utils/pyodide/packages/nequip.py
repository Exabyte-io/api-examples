import os
import sys
import types
from importlib import import_module

import torch

from ...primitive.environment import is_pyodide_environment
from .torch import _ensure_omegaconf_stub, _make_stub_module, _matscipy_neighbour_list_compat


def apply_patches():
    patch_nequip_deps()


def patch_nequip_deps():
    """
    Stub heavy dependencies required by NequIP but not needed for inference.

    Stubs: hydra, lightning, pytorch_lightning, torchmetrics, lmdb, matscipy.
    Patches e3nn and torch.jit for Pyodide compatibility.
    Sets NEQUIP_NL=ase to use ASE neighbor lists instead of matscipy.
    """
    os.environ["NEQUIP_NL"] = "ase"
    _patch_lightning_utilities()
    _patch_lightning()
    _patch_pytorch_lightning()
    _patch_torchmetrics()
    _patch_packaging()
    _patch_lmdb()
    _patch_hydra()
    _ensure_omegaconf_stub()
    _patch_matscipy()
    _patch_e3nn()
    _patch_torch_jit()
    _patch_e3nn_jit()
    _patch_tqdm()

    print("✓ NequIP dependency stubs applied")


def _patch_lightning_utilities():
    _make_stub_module("lightning_utilities", submodules=["core", "core.rank_zero"])

    def _rank_prefixed_message(msg, rank=None):
        return msg

    def _rank_zero_only(fn):
        return fn

    rank_zero_mod = sys.modules["lightning_utilities.core.rank_zero"]
    rank_zero_mod.rank_prefixed_message = _rank_prefixed_message
    rank_zero_mod.rank_zero_only = _rank_zero_only


def _patch_lightning():
    _make_stub_module(
        "lightning",
        submodules=[
            "pytorch",
            "pytorch.utilities",
            "pytorch.utilities.seed",
            "pytorch.utilities.warnings",
            "pytorch.callbacks",
        ],
    )

    class _IsolateRng:
        def __enter__(self):
            return self

        def __exit__(self, *a):
            pass

    def _seed_everything(seed=None, workers=False, **kwargs):
        pass

    sys.modules["lightning.pytorch.utilities.seed"].isolate_rng = lambda: _IsolateRng()
    sys.modules["lightning.pytorch"].seed_everything = _seed_everything

    class _PossibleUserWarning(UserWarning):
        pass

    sys.modules["lightning.pytorch.utilities.warnings"].PossibleUserWarning = _PossibleUserWarning

    class _LightningModule(torch.nn.Module):
        def __init__(self, *a, **k):
            super().__init__()

        def log(self, *a, **k):
            pass

    sys.modules["lightning.pytorch"].LightningModule = _LightningModule
    sys.modules["lightning"].pytorch = sys.modules["lightning.pytorch"]
    sys.modules["lightning.pytorch.callbacks"].Callback = type("Callback", (), {})


def _patch_pytorch_lightning():
    _make_stub_module("pytorch_lightning", submodules=["utilities", "utilities.seed"])

    class _IsolateRng:
        def __enter__(self):
            return self

        def __exit__(self, *a):
            pass

    sys.modules["pytorch_lightning.utilities.seed"].isolate_rng = lambda: _IsolateRng()


def _patch_torchmetrics():
    torchmetrics_mod = _make_stub_module("torchmetrics")

    class _Metric(torch.nn.Module):
        def __init__(self, *a, **k):
            super().__init__()

        def add_state(self, name, default=None, dist_reduce_fx=None, **k):
            if default is not None:
                setattr(self, name, default)

    torchmetrics_mod.Metric = _Metric
    torchmetrics_mod.MeanMetric = _Metric


def _patch_packaging():
    if "packaging" in sys.modules:
        return
    try:
        import packaging.version  # noqa: F401
    except (ImportError, ModuleNotFoundError):
        _make_stub_module("packaging", submodules=["version"])

        class _Version:
            def __init__(self, v):
                self._v = str(v)
                parts = self._v.split(".")
                self.major = int(parts[0]) if parts else 0
                self.minor = int(parts[1]) if len(parts) > 1 else 0

            def __lt__(self, other):
                return (self.major, self.minor) < (other.major, other.minor)

            def __ge__(self, other):
                return not self.__lt__(other)

            def __repr__(self):
                return f"Version('{self._v}')"

        sys.modules["packaging.version"].Version = _Version
        sys.modules["packaging.version"].parse = lambda v: _Version(v)


def _patch_lmdb():
    if "lmdb" not in sys.modules:
        sys.modules["lmdb"] = types.ModuleType("lmdb")


def _patch_hydra():
    if "hydra" in sys.modules:
        return
    _make_stub_module(
        "hydra",
        submodules=[
            "core",
            "core.global_hydra",
            "utils",
            "_internal",
            "_internal.instantiate",
            "_internal.instantiate._instantiate2",
        ],
    )

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
    sys.modules["hydra._internal.instantiate._instantiate2"].InstantiationException = type(
        "InstantiationException", (Exception,), {}
    )

    def _hydra_get_target(target_str):
        mod_path, name = target_str.rsplit(".", 1)
        return getattr(import_module(mod_path), name)

    sys.modules["hydra.utils"].get_method = _hydra_get_target
    sys.modules["hydra.utils"].get_class = _hydra_get_target


def _patch_matscipy():
    if "matscipy" in sys.modules:
        return
    matscipy_mod = types.ModuleType("matscipy")
    matscipy_mod.__path__ = []
    matscipy_mod.__package__ = "matscipy"
    matscipy_neighbours = types.ModuleType("matscipy.neighbours")
    matscipy_neighbours.neighbour_list = _matscipy_neighbour_list_compat
    matscipy_mod.neighbours = matscipy_neighbours
    sys.modules["matscipy"] = matscipy_mod
    sys.modules["matscipy.neighbours"] = matscipy_neighbours


def _patch_e3nn():
    try:
        import e3nn

        e3nn._SO3_INITIALIZED = True
        if hasattr(e3nn, "_OPT_DEFAULTS"):
            e3nn._OPT_DEFAULTS["jit_mode"] = "eager"
    except Exception:
        pass


def _patch_torch_jit():
    if not hasattr(torch.jit, "_original_script"):

        def _noop_script(obj=None, *a, **k):
            if obj is not None:
                return obj
            return lambda fn: fn

        torch.jit.script = _noop_script


def _patch_e3nn_jit():
    try:
        import e3nn.util.jit

        e3nn.util.jit.compile_mode = lambda mode: lambda cls: cls
    except Exception:
        pass


def _patch_tqdm():
    if "tqdm" in sys.modules:
        return
    tqdm_mod = _make_stub_module("tqdm", submodules=["auto", "std"])

    def _tqdm_passthrough(iterable=None, *a, **k):
        return iterable if iterable is not None else iter([])

    tqdm_mod.tqdm = _tqdm_passthrough
    sys.modules["tqdm.auto"].tqdm = _tqdm_passthrough
    sys.modules["tqdm.std"].tqdm = _tqdm_passthrough


MODEL_PATHS_MAP = {
    "oam_s": "/drive/packages/models/nequip-oam-s-config-sd.pth",
}


def get_nequip_model_pyodide(model: str, device="cpu"):
    if model not in MODEL_PATHS_MAP:
        raise ValueError(f"Invalid model name: {model}. Valid options are: {list(MODEL_PATHS_MAP.keys())}")

    nequip_model = load_nequip_model(MODEL_PATHS_MAP[model])
    return create_nequip_calculator_from_model(nequip_model, device=device)


def create_nequip_calculator_from_model(nequip_model, device="cpu"):
    r_max = float(nequip_model.metadata["r_max"])
    type_names = nequip_model.metadata["type_names"].split(" ")
    nequip_ase = import_module("nequip.ase")
    transforms = import_module("nequip.data.transforms")

    return nequip_ase.NequIPCalculator(
        model=nequip_model,
        device=device,
        transforms=[
            transforms.ChemicalSpeciesToAtomTypeMapper(type_names),
            transforms.NeighborListTransform(r_max=r_max),
        ],
    )


def create_nequip_calculator(model="oam_s", device="cpu", model_path=None, checkpoint=None):
    if is_pyodide_environment():
        return get_nequip_model_pyodide(model=model, device=device)

    resolved_model_path = model_path or checkpoint

    nequip_model = load_nequip_model(str(resolved_model_path))
    return create_nequip_calculator_from_model(nequip_model, device=device)


def load_nequip_model(checkpoint_path):
    data = torch.load(checkpoint_path, map_location="cpu", weights_only=False)

    import nequip.utils.global_state as global_state

    def _pyodide_set_global_state(allow_tf32=False, warn_on_override=False):
        if not global_state._GLOBAL_STATE_INITIALIZED:
            torch.set_default_dtype(torch.float64)
            try:
                import e3nn

                e3nn.set_optimization_defaults(
                    specialized_code=True,
                    optimize_einsums=True,
                    jit_script_fx=False,
                )
            except Exception:
                pass
            global_state._GLOBAL_STATE_INITIALIZED = True
        global_state._latest_global_config["allow_tf32"] = allow_tf32

    global_state.set_global_state = _pyodide_set_global_state
    global_state.set_global_state(allow_tf32=False)

    from nequip.model import FullNequIPGNNModel

    model = FullNequIPGNNModel(
        seed=0,
        model_dtype="float32",
        r_max=data["r_max"],
        type_names=data["type_names"],
        irreps_edge_sh=data["irreps_edge_sh"],
        type_embed_num_features=data["type_embed_num_features"],
        feature_irreps_hidden=data["feature_irreps_hidden"],
        radial_mlp_depth=data["radial_mlp_depth"],
        radial_mlp_width=data["radial_mlp_width"],
        avg_num_neighbors=data["avg_num_neighbors"],
        per_type_energy_scales=data["per_type_energy_scales"],
        per_type_energy_shifts=data["per_type_energy_shifts"],
        polynomial_cutoff_p=data.get("polynomial_cutoff_p", 6),
    )

    if data.get("has_zbl", False):
        from nequip.nn.pair_potential import ZBL

        seq_net = model.model.func
        zbl = ZBL(
            type_names=data["type_names"],
            chemical_species=data["type_names"],
            units="metal",
            irreps_in=seq_net.irreps_out,
        )
        seq_net.insert(name="pair_potential", module=zbl, before="total_energy_sum")

    model.load_state_dict(data["state_dict"], strict=True)
    model.eval()

    print(f"✓ NequIP model loaded ({sum(p.numel() for p in model.parameters()):,} parameters)")
    return model
