from .generic import ensure_packages_installed, update_json_file_kwargs

def setup_colab_runtime_environment(environment_variables_config):
    """
    This function setups up the runtime environment for running the exabyte-api-examples
    notebooks within Google Colaboratory

    Args:
        environment_variables_config (dict): Contains the variables needed to properly setup the
            a runtime_environment

    Returns
        None
    """
    ensure_packages_installed(environment_variables_config)
    kwargs = {key: environment_variables_config[key] for key in ["ACCOUNT_ID", "AUTH_TOKEN", "MATERIALS_PROJECT_API_KEY", "ORGANIZATION_ID"]}
    update_json_file_kwargs(**kwargs)
