import os

from .enums import EnvironmentsEnum

# default value for env.HOME from https://pyodide.org/en/stable/usage/api/js-api.html
ENVIRONMENT = EnvironmentsEnum.PYODIDE if os.environ.get("HOME") == "/home/pyodide" else EnvironmentsEnum.PYTHON


def is_pyodide_environment() -> bool:
    return ENVIRONMENT == EnvironmentsEnum.PYODIDE
