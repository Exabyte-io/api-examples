import inspect
import os

from mat3ra.api_client import ACCESS_TOKEN_ENV_VAR

from .environment import ENVIRONMENT, EnvironmentsEnum
from .io import get_data
from .jupyterlite.api.auth import authenticate_jupyterlite, authenticate_oidc

REFRESH_TOKEN_ENV_VAR = "OIDC_REFRESH_TOKEN"


async def authenticate(force=False, globals_dict=None):
    if globals_dict is None:
        frame = inspect.currentframe()
        try:
            globals_dict = frame.f_back.f_globals  # type: ignore
        finally:
            del frame
    if ENVIRONMENT == EnvironmentsEnum.PYODIDE:
        get_data("data_from_host", globals_dict)
    data_from_host = globals_dict.get("data_from_host")
    if data_from_host:
        await authenticate_jupyterlite(data_from_host)
    elif ACCESS_TOKEN_ENV_VAR not in os.environ or force:
        await authenticate_oidc()
