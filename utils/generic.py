# This module defines a set of common functions which are used in other examples.
import json
import os
import uuid
from types import SimpleNamespace
from typing import Union

from IPython.display import HTML, display
from pandas import DataFrame
from pandas.io.formats.style import Styler

from . import settings


def update_json_file_kwargs(path_to_json_file: str = "settings.json", **kwargs) -> None:
    """
    This function updates settings.json for a given kwargs if kwargs
    contains variables different from those already in json

    Args:
        path_to_json_file (str): the path to the json file to be updated
        **kwargs (dict): A dict of keyword arguments

    Returns:
        None
    """

    # 1. Assert the json file is where we think it is
    assert os.path.isfile(path_to_json_file)

    # 2. Load settings.json
    with open(path_to_json_file) as settings_json_file:
        variables = json.load(settings_json_file)

    # 3. Update json file if kwargs contains new variables
    if kwargs != variables:
        updated_variables = {**variables, **kwargs}
        with open(path_to_json_file, "w") as settings_json_file:
            json.dump(updated_variables, settings_json_file, indent=4)


# DISPLAY UTILITIES


def dataframe_to_html(df: DataFrame, text_align: str = "center") -> Styler:
    """
    Converts Pandas dataframe to HTML.
    See https://pandas.pydata.org/pandas-docs/stable/style.html for more information about styling.

    Args:
        df (pd.DataFrame): Pandas dataframe.
        text_align (str): text align. Defaults to center.
    """
    styles = [
        dict(selector="th", props=[("text-align", text_align)]),
        dict(selector="td", props=[("text-align", text_align)]),
    ]
    return df.style.set_table_styles(styles)


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


# Helper function to convert dictionaries to SimpleNamespace objects for dot notation access
def dict_to_namespace(obj):
    if isinstance(obj, dict):
        return SimpleNamespace(**{k: dict_to_namespace(v) for k, v in obj.items()})
    elif isinstance(obj, list):
        return [dict_to_namespace(item) for item in obj]
    else:
        return obj
