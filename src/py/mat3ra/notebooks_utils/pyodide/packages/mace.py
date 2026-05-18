from mace.calculators import MACECalculator

MODEL_PATHS_MAP = {
    "small": "/drive/packages/models/2023-12-10-mace-128-L0_energy_epoch-249.model",
    "medium": "/drive/packages/models/2023-12-03-mace-128-L1_epoch-199.model",
    "large": "/drive/packages/models/MACE_MPtrj_2022.9.model",
}


def get_mace_model_pyodide(model_name):
    if model_name not in MODEL_PATHS_MAP:
        raise ValueError(f"Invalid model name: {model_name}. Valid options are: {list(MODEL_PATHS_MAP.keys())}")
    model_path = MODEL_PATHS_MAP[model_name]
    return MACECalculator(model_path=model_path, dispersion=False, default_dtype="float32", device="cpu")
