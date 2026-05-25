from mattersim.forcefield import MatterSimCalculator

from ...primitive.environment import is_pyodide_environment

MODEL_PATHS_MAP = {
    "1m": "/drive/packages/models/mattersim-v1.0.0-1M.pth",
}


def get_mattersim_model_pyodide(model: str, device="cpu", **kwargs):
    if model not in MODEL_PATHS_MAP:
        raise ValueError(f"Invalid model name: {model}. Valid options are: {list(MODEL_PATHS_MAP.keys())}")
    return MatterSimCalculator.from_checkpoint(load_path=MODEL_PATHS_MAP[model], device=device, **kwargs)


def create_mattersim_calculator(model="1m", device="cpu", model_path=None, checkpoint=None, **kwargs):
    if is_pyodide_environment():
        return get_mattersim_model_pyodide(model=model, device=device, **kwargs)

    resolved_model_path = model_path or checkpoint

    return MatterSimCalculator.from_checkpoint(load_path=str(resolved_model_path), device=device, **kwargs)
