import json
import os
import uuid
from typing import Union

from IPython.display import HTML, display

from . import settings


def display_JSON(
    obj: Union[dict, list], interactive_viewer: bool = settings.use_interactive_JSON_viewer, level: int = 2
) -> None:
    """
    Displays JSON, either interactively or via a text dump to Stdout.

    The interactive viewer is based on https://github.com/mljar/mercury/blob/main/mercury/widgets/json.py.

    Args:
        obj (dict): Object to display as nicely-formatted JSON
        interactive_viewer (bool): Whether to use the interactive viewer or not
        level (int): The level to which the JSON should be expanded by default
    """
    if interactive_viewer:
        if isinstance(obj, (dict, list)):
            json_str = json.dumps(obj)
        else:
            json_str = obj

        id = str(uuid.uuid4())

        web_dir = os.path.join(os.path.dirname(__file__), "web")

        with open(os.path.join(web_dir, "renderjson.css")) as fp:
            css = fp.read()

        with open(os.path.join(web_dir, "renderjson.js")) as fp:
            js = fp.read()

        display(HTML(f'<style>{css}</style><div id="{id}"></div>'))
        display(
            HTML(
                f"<script>{js} "
                f"renderjson.set_show_to_level({str(level)}); "
                f'renderjson.set_icons("▸","▾"); '
                f'document.getElementById("{id}").appendChild(renderjson({json_str}))</script>'
            )
        )
    else:
        print(json.dumps(obj, indent=4))
