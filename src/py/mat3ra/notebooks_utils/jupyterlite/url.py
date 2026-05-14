from typing import Union

from ..enums import EnvironmentsEnum
from ..environment import ENVIRONMENT


def read_from_url_python(url: str, as_bytes: bool = False) -> Union[str, bytes]:
    """
    Fetch, read, and decode content from a URL in a Python environment.

    Args:
        url (str): The URL to fetch from.
        as_bytes (bool): Whether to return the content as bytes.

    Returns:
        str: The content as a string or bytes.
    """
    import urllib.request

    with urllib.request.urlopen(url) as response:
        body = response.read()
        if as_bytes:
            return body
        return body.decode("utf-8")


async def read_from_url_pyodide(url: str, as_bytes: bool = False) -> Union[str, bytes]:
    """
    Fetch and read content from a URL in a Pyodide environment.

    Args:
        url (str): The URL to fetch from.
        as_bytes (bool): Whether to return the content as bytes.

    Returns:
        str: The content as a string or bytes.
    """
    # `http` is a Pyodide module that will be installed in the Pyodide environment by default.
    from pyodide.http import pyfetch  # type: ignore

    # Per https://developer.mozilla.org/en-US/docs/Web/API/Fetch_API/Using_Fetch
    response = await pyfetch(url)
    if as_bytes:
        return await response.bytes()
    return await response.string()


async def read_from_url(url: str, as_bytes: bool = False) -> Union[str, bytes]:
    """
    Read content from a URL, handling both Python and Pyodide environments.

    Args:
        url (str): The URL to fetch from.
        as_bytes (bool): Whether to return the content as

    Returns:
        str: The content as a string or bytes.
    """
    response = None
    if ENVIRONMENT == EnvironmentsEnum.PYODIDE:
        response = await read_from_url_pyodide(url, as_bytes)
    elif ENVIRONMENT == EnvironmentsEnum.PYTHON:
        response = read_from_url_python(url, as_bytes)

    return response
