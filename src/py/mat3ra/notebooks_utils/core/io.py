import json
import os
import urllib.request
from typing import Any, Dict, Optional

from ..primitive.logger import log
from ..settings import UPLOADS_FOLDER


def get_data_python(key: str, globals_dict: Optional[Dict] = None):
    """
    Read data from the `uploads` folder in a JupyterLab environment.

    Args:
        key (str): The name under which data is expected to be received.
        globals_dict (dict, optional): A dictionary to store the received data. Defaults to None.
    """
    try:
        data_from_host = []
        index = 0
        for filename in sorted(os.listdir(UPLOADS_FOLDER)):
            if filename.endswith(".json"):
                with open(os.path.join(UPLOADS_FOLDER, filename), "r") as file:
                    data = json.load(file)
                name = os.path.splitext(filename)[0]
                log(f"{index}: {name}")
                index += 1
                data_from_host.append(data)
        if globals_dict is not None:
            globals_dict[key] = data_from_host
        return data_from_host
    except FileNotFoundError:
        print("No data found in the 'uploads' folder.")


def set_data_python(key: str, value: Any):
    """
    Write data to the `uploads` folder in a JupyterLab environment.

    Args:
        key (str): The name under which data will be written.
        value (Any): The value to write to the `uploads` folder.
    """
    if not os.path.exists(UPLOADS_FOLDER):
        os.makedirs(UPLOADS_FOLDER)
    for item in value:
        safe_name = item["name"].replace("%", "pct").replace("/", ":")
        file_path = os.path.join(UPLOADS_FOLDER, f"{safe_name}.json")
        with open(file_path, "w") as file:
            json.dump(item, file)
        log(f"Data for {key} written to {file_path}")


def read_from_url_python(url: str, as_bytes: bool = False):
    """
    Fetch, read, and decode content from a URL in a Python environment.

    Args:
        url (str): The URL to fetch from.
        as_bytes (bool): Whether to return the content as bytes.

    Returns:
        str or bytes: The content.
    """
    with urllib.request.urlopen(url) as response:
        body = response.read()
        if as_bytes:
            return body
        return body.decode("utf-8")
