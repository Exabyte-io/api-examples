import json
import os
from types import SimpleNamespace


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


# Helper function to convert dictionaries to SimpleNamespace objects for dot notation access
def dict_to_namespace(obj):
    if isinstance(obj, dict):
        return SimpleNamespace(**{k: dict_to_namespace(v) for k, v in obj.items()})
    elif isinstance(obj, list):
        return [dict_to_namespace(item) for item in obj]
    else:
        return obj
