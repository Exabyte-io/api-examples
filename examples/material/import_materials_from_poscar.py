#!/usr/bin/env python
# coding: utf-8

# <a href="https://colab.research.google.com/github/Exabyte-io/api-examples/blob/dev/examples/material/import_materials_from_poscar.ipynb" target="_parent">
# <img alt="Open in Google Colab" src="https://user-images.githubusercontent.com/20477508/128780728-491fea90-9b23-495f-a091-11681150db37.jpeg" width="150" border="0">
# </a>

# # Overview
# 
# This example demonstrates how to import a material from a POSCAR file via [Material](https://docs.mat3ra.com/api/Material/post_materials_import) endpoints.

# # Complete Authorization Form and Initialize Settings
# 
# This will also determine environment and set all environment variables. We determine if we are using Jupyter Notebooks or Google Colab to run this tutorial.
# 
# If you are running this notebook from Google Colab, Colab takes ~1 min to execute the following cell.
# 
# ACCOUNT_ID and AUTH_TOKEN - Authentication parameters needed for when making requests to [Mat3ra.com's API Endpoints](https://docs.mat3ra.com/rest-api/endpoints/).
# 
# MATERIALS_PROJECT_API_KEY - Authentication parameter needed for when making requests to [Material Project's API](https://materialsproject.org/open)
# 
# ORGANIZATION_ID - Authentication parameter needed for when working with collaborative accounts https://docs.mat3ra.com/collaboration/organizations/overview/
# 
# > <span style="color: orange">**NOTE**</span>: If you are running this notebook from Jupyter, the variables ACCOUNT_ID, AUTH_TOKEN, MATERIALS_PROJECT_API_KEY, and ORGANIZATION_ID should be set in the file [settings.json](../../utils/settings.json) if you need to use these variables. To obtain API token parameters, please see the following link to the documentation explaining how to get them: https://docs.mat3ra.com/accounts/ui/preferences/api/

# In[ ]:


# @title Authorization Form
ACCOUNT_ID = "ACCOUNT_ID"  # @param {type:"string"}
AUTH_TOKEN = "AUTH_TOKEN"  # @param {type:"string"}
MATERIALS_PROJECT_API_KEY = "MATERIALS_PROJECT_API_KEY"  # @param {type:"string"}
ORGANIZATION_ID = "ORGANIZATION_ID"  # @param {type:"string"}

import os

if "COLAB_JUPYTER_IP" in os.environ:
    os.environ.update(
        dict(
            ACCOUNT_ID=ACCOUNT_ID,
            AUTH_TOKEN=AUTH_TOKEN,
            MATERIALS_PROJECT_API_KEY=MATERIALS_PROJECT_API_KEY,
            ORGANIZATION_ID=ORGANIZATION_ID,
        )
    )

    get_ipython().system('GIT_BRANCH="dev"; export GIT_BRANCH; export IS_USING_GIT_LFS=true; curl -s "https://raw.githubusercontent.com/Exabyte-io/api-examples/${GIT_BRANCH}/scripts/env.sh" | bash')
    from examples.utils.notebook import get_notebook_info

    os.chdir(os.path.join("api-examples", os.path.dirname(get_notebook_info()["notebook_path"])))


# # Imports

# In[ ]:


from utils.settings import ENDPOINT_ARGS
from utils.generic import display_JSON

from exabyte_api_client.endpoints.materials import MaterialEndpoints


# ## Set Parameters
# 
# - **NAME**: material name
# - **POSCAR_PATH**: absolute path to the POSCAR file

# In[ ]:


NAME = "My Material"
POSCAR_PATH = "../assets/mp-978534.poscar"


# ## Import material
# 
# Initialize `MaterialEndpoints` class and call `import_from_file` function to import the material.

# In[ ]:


content = ""
with open(POSCAR_PATH) as f:
    content = f.read()
    print(content)

endpoint = MaterialEndpoints(*ENDPOINT_ARGS)
material = endpoint.import_from_file(NAME, content)


# ## Print imported material
# 
# Print the list of imported materials in pretty JSON below.

# In[ ]:


display_JSON(material)

