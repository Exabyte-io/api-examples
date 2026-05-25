from mace.calculators import MACECalculator, mace_mp

from ...primitive.environment import is_pyodide_environment

MODEL_PATHS_MAP = {
    "small": "/drive/packages/models/2023-12-10-mace-128-L0_energy_epoch-249.model",
    "medium": "/drive/packages/models/2023-12-03-mace-128-L1_epoch-199.model",
    "large": "/drive/packages/models/MACE_MPtrj_2022.9.model",
}


def get_mace_model_pyodide(model: str, dispersion=False, default_dtype="float32", device="cpu", **kwargs):
    if model not in MODEL_PATHS_MAP:
        raise ValueError(f"Invalid model name: {model}. Valid options are: {list(MODEL_PATHS_MAP.keys())}")
    model_path = MODEL_PATHS_MAP[model]
    return MACECalculator(
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

    return mace_mp(
        model=model,
        dispersion=dispersion,
        default_dtype=default_dtype,
        device=device,
        **kwargs,
    )
