# The variables defined below in environment_variables_config are defined in the first code cell
# in either the Google Colab or Jupyter notebook.
#
# In Google Colab, these variables must be filled out in the 'Authorization Form' section of the notebook.
#
# If using Jupyter, these variables can be left to their default values in the code cell, but the user
# should change these values in the settings.json file located in the examples folder.

import os
import re
import sys
from urllib.parse import unquote, urlparse

import requests
from IPython import get_ipython

environment_variables_config = {
    "ACCOUNT_ID": os.getenv("ACCOUNT_ID", ""),
    "AUTH_TOKEN": os.getenv("AUTH_TOKEN", ""),
    "MATERIALS_PROJECT_API_KEY": os.getenv("MATERIALS_PROJECT_API_KEY", ""),
    "ORGANIZATION_ID": os.getenv("ORGANIZATION_ID", ""),
}


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
    notebook_environment = environment_variables_config["notebook_environment"]
    if notebook_environment == "Colab":
        from utils.colab import setup_colab_runtime_environment

        setup_colab_runtime_environment(environment_variables_config)
    else:
        from utils.generic import ensure_packages_are_installed

        ensure_packages_are_installed(notebook_environment)


def get_notebook_info():
    """
    Get the name of a currently running notebook in Google Colab.
    Args:
        None
    Return:
        filename
    """
    # ip = socket.gethostbyname(socket.gethostname())  # 172.28.0.12
    ip = os.getenv("COLAB_JUPYTER_IP")
    response = requests.get(f"http://{ip}:9000/api/sessions").json()[0]

    notebook_name = response["name"]
    path = urlparse(unquote(response["path"]).split("=")[1]).path
    parsed = re.findall("(.*)/blob/(.*)/examples/(.*)", path)[0]
    github_org_repo = parsed[0][1:]  # remove leading /
    branch_name = parsed[1]
    notebook_path = f"examples/{parsed[2]}"

    print(notebook_path)

    return dict(
        notebook_name=notebook_name,
        notebook_path=notebook_path,
        branch_name=branch_name,
        github_org_repo=github_org_repo,
    )


def execute():
    """
    Main execution function. This function determines and sets the runtime environment
    for a given notebook frontend, such as Jupyter Notebooks or Google Colab.
    Args:
        None
    Return:
        None
    """

    if "is_setup_executed" not in os.environ:
        ip = get_ipython()
        if "google.colab" in str(ip):
            environment_variables_config.update({"notebook_environment": "Colab"})

            notebook_info = get_notebook_info()
            notebook_path = notebook_info["notebook_path"]
            github_org_repo = notebook_info["github_org_repo"]
            repo_name = github_org_repo.split("/")[1]

            # Go to the folder in the repo where one would be if this was in local Jupyter.
            # The repo should be clonned via the notebook cell earlier.
            os.chdir(os.path.join(repo_name, os.path.dirname(notebook_path)))

        elif "ZMQInteractiveShell" in str(ip):
            environment_variables_config.update({"notebook_environment": "Jupyter"})
        else:
            environment_variables_config.update({"notebook_environment": ""})

        module_path = os.path.abspath(os.path.join(".."))
        if module_path not in sys.path:
            sys.path.append(module_path)
        set_notebook_environment(environment_variables_config)


execute()
os.environ.update({"is_setup_executed": "True"})
