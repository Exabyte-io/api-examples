from fairchem.core import FAIRChemCalculator
from fairchem.core.units.mlip_unit import load_predict_unit

from ...primitive.environment import is_pyodide_environment

MODEL_PATHS_MAP = {
    "f16": "/drive/packages/models/uma-s-1p1-f16.pt",
    "int8": "/drive/packages/models/uma-s-1p1-int8.pt",
}


def get_uma_model_pyodide(model: str, task_name="omat", device="cpu", **kwargs):
    if model not in MODEL_PATHS_MAP:
        raise ValueError(f"Invalid model name: {model}. Valid options are: {list(MODEL_PATHS_MAP.keys())}")
    predictor = load_predict_unit(MODEL_PATHS_MAP[model], device=device, **kwargs)
    return FAIRChemCalculator(predictor, task_name=task_name)


def create_uma_calculator(model="f16", task_name="omat", device="cpu", model_path=None, checkpoint=None, **kwargs):
    if is_pyodide_environment():
        return get_uma_model_pyodide(model=model, task_name=task_name, device=device, **kwargs)

    resolved_model_path = model_path or checkpoint

    predictor = load_predict_unit(str(resolved_model_path), device=device, **kwargs)
    return FAIRChemCalculator(predictor, task_name=task_name)
