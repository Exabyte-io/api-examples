from types import SimpleNamespace


# TODO: dict_to_namespace_recursive move to mat3ra.utils
# Helper function to convert dictionaries to SimpleNamespace objects for dot notation access
def dict_to_namespace(obj):
    if isinstance(obj, dict):
        return SimpleNamespace(**{k: dict_to_namespace(v) for k, v in obj.items()})
    elif isinstance(obj, list):
        return [dict_to_namespace(item) for item in obj]
    else:
        return obj
