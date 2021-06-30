environment_variables_config = {'ACCOUNT_ID': ACCOUNT_ID,
                                'AUTH_TOKEN': AUTH_TOKEN,
                                'MATERIALS_PROJECT_API_KEY': MATERIALS_PROJECT_API_KEY,
                                'ORGANIZATION_ID': ORGANIZATION_ID}

def set_notebook_environment(environment_variables_config):
    """
    This function sets the notebook environment by calling to install the needed packages
    and setting the variables in settings.json (if needed)
    Args:
        environment_variables_config (dict): contains key value pairs needed to set up a
            certain notebook (or frontend) runtime.
            Ex) environment_variables_config['ACCOUNT_ID'] = ACCOUNT_ID
                environment_variables_config[notebook_environment] = "Jupyter"

    Return:
        None
    """
    notebook_environment = environment_variables_config['notebook_environment']
    if notebook_environment == 'Colab':
        from utils.colab import setup_colab_runtime_environment
        setup_colab_runtime_environment(environment_variables_config)
    else:
        from utils.generic import ensure_packages_are_installed
        ensure_packages_are_installed(notebook_environment)


def execute():
    """
    Main execution function. This function determines and sets the runtime environment
    for a given notebook frontend, such as Jupyter Notebooks or Google Colab.
    Args:
        None
    Return:
        None
    """

    if 'google.colab' in str(get_ipython()):
        environment_variables_config.update({'notebook_environment': 'Colab'})
        get_ipython().system('git clone -b feature/SOF-4400-skinny-req https://github.com/Exabyte-io/exabyte-api-examples.git')
        from google.colab import _message
        notebook_name = _message.blocking_request('get_ipynb')['ipynb']['metadata']['colab']['name']
        notebook_path = glob.glob('**/'+notebook_name, recursive=True)[0][0:-len(notebook_name)]
        os.chdir(notebook_path) # go to the folder in the repo where one would be if this was in local Jupyter
    elif 'ZMQInteractiveShell' in str(get_ipython()):
        environment_variables_config.update({'notebook_environment': 'Jupyter'})
    else:
        environment_variables_config.update({'notebook_environment': ''})

    module_path = os.path.abspath(os.path.join('..'))
    if module_path not in sys.path: sys.path.append(module_path)
    set_notebook_environment(environment_variables_config)

execute()
