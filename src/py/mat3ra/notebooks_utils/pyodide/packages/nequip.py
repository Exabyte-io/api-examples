from nequip.ase import NequIPCalculator
from nequip.data.transforms import ChemicalSpeciesToAtomTypeMapper, NeighborListTransform

from ...primitive.environment import is_pyodide_environment
from .torch import load_nequip_model

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

    return NequIPCalculator(
        model=nequip_model,
        device=device,
        transforms=[
            ChemicalSpeciesToAtomTypeMapper(type_names),
            NeighborListTransform(r_max=r_max),
        ],
    )


def create_nequip_calculator(model="oam_s", device="cpu", model_path=None, checkpoint=None):
    if is_pyodide_environment():
        return get_nequip_model_pyodide(model=model, device=device)

    resolved_model_path = model_path or checkpoint

    nequip_model = load_nequip_model(str(resolved_model_path))
    return create_nequip_calculator_from_model(nequip_model, device=device)
