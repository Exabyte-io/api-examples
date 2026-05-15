from enum import Enum


class EnvironmentsEnum(Enum):
    PYODIDE = "pyodide"
    PYTHON = "python"


class SeverityLevelEnum(Enum):
    INFO = "INFO"
    WARNING = "WARNING"
    ERROR = "ERROR"
