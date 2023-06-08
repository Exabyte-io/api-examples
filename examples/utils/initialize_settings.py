# The variables defined below in environment_variables_config are defined in the first code cell
# in either the Google Colab or Jupyter notebook.
#
# In Google Colab, these variables must be filled out in the 'Authorization Form' section of the notebook.
#
# If using Jupyter, these variables can be left to their default values in the code cell, but the user
# should change these values in the settings.json file located in the examples folder.

import glob
import os
import socket
import sys

import requests
from IPython import get_ipython

environment_variables_config = {'ACCOUNT_ID': ACCOUNT_ID,  # noqa F821
                                'AUTH_TOKEN': AUTH_TOKEN,  # noqa F821
                                'MATERIALS_PROJECT_API_KEY': MATERIALS_PROJECT_API_KEY,  # noqa F821
                                'ORGANIZATION_ID': ORGANIZATION_ID}  # noqa F821


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


def get_notebook_name():
    """
    Get the name of a currently running notebook in Google Colab.
    Args:
        None
    Return:
        filename
    """
    ip = socket.gethostbyname(socket.gethostname())  # 172.28.0.12
    filename = requests.get(f"http://{ip}:9000/api/sessions").json()[0]["name"]
    return filename


def execute():
    """
    Main execution function. This function determines and sets the runtime environment
    for a given notebook frontend, such as Jupyter Notebooks or Google Colab.
    Args:
        None
    Return:
        None
    """

    if 'is_setup_executed' not in os.environ:
        ip = get_ipython()
        if 'google.colab' in str(ip):
            branch = os.getenv("GIT_BRANCH", "dev")  # a way to inject a branch of interest from Colab if run via a PR.
            environment_variables_config.update({'notebook_environment': 'Colab'})
            ip.system(f"git clone --single-branch -b {branch} https://github.com/Exabyte-io/exabyte-api-examples.git")
            notebook_name = get_notebook_name()
            notebook_path = glob.glob(f"**/{notebook_name}", recursive=True)[0][0:-len(notebook_name)]
            os.chdir(notebook_path)  # go to the folder in the repo where one would be if this was in local Jupyter
        elif 'ZMQInteractiveShell' in str(ip):
            environment_variables_config.update({'notebook_environment': 'Jupyter'})
        else:
            environment_variables_config.update({'notebook_environment': ''})

        module_path = os.path.abspath(os.path.join('..'))
        if module_path not in sys.path:
            sys.path.append(module_path)
        set_notebook_environment(environment_variables_config)


execute()
os.environ.update({'is_setup_executed': 'True'})
