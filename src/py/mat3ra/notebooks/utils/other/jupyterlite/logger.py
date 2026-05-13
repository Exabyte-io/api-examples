import inspect
import os
from typing import Optional

from .enums import SeverityLevelEnum


def log(message: str, level: Optional[SeverityLevelEnum] = None, force_verbose: Optional[bool] = None):
    """
    Log a message based on the VERBOSE flag in the caller's globals().

    Args:
        message (str): The message to log.
        level (SeverityLevelEnum, optional): The severity level of the message (e.g., INFO, WARNING, ERROR).
        force_verbose (bool, optional): If True, log the message regardless of the VERBOSE flag in globals().
    """
    if force_verbose is True:
        should_log = True
    elif force_verbose is False:
        should_log = False
    else:
        # Inspect the caller's globals to get VERBOSE flag
        frame = inspect.currentframe()
        try:
            caller_frame = frame.f_back  # type: ignore
            caller_globals = caller_frame.f_globals  # type: ignore
            should_log = caller_globals.get("VERBOSE", os.environ.get("VERBOSE", True))
        finally:
            del frame  # Avoid reference cycles
    if should_log:
        if level is None:
            print(message)
        else:
            print(f"{level.value}: {message}")
