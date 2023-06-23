#!/usr/bin/env python
# coding: utf-8

# <a href="https://colab.research.google.com/github/Exabyte-io/api-examples/blob/dev/examples/material/create_material.ipynb" target="_parent">
# <img alt="Open in Google Colab" src="https://user-images.githubusercontent.com/20477508/128780728-491fea90-9b23-495f-a091-11681150db37.jpeg" width="150" border="0">
# </a>

# # Overview
# 
# In this example we create a material from a JSON config with [tags](https://docs.mat3ra.com/entities-general/data/#tags) to identify the material.

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

    get_ipython().system('GIT_BRANCH="dev"; export GIT_BRANCH; curl -s "https://raw.githubusercontent.com/Exabyte-io/api-examples/${GIT_BRANCH}/scripts/env.sh" | bash')


# # Imports

# In[ ]:


from utils.settings import ENDPOINT_ARGS
from utils.generic import display_JSON

from exabyte_api_client.endpoints.materials import MaterialEndpoints


# ## Create material config
# 
# Create material config in JSON format. See [Material](https://docs.mat3ra.com/api/Material/put_materials_create) endpoint for more information about material config format.

# In[ ]:


CONFIG = {
    "name": "TEST MATERIAL",
    "basis": {
        "elements": [{"id": 1, "value": "Si"}, {"id": 2, "value": "Si"}],
        "coordinates": [{"id": 1, "value": [0, 0, 0]}, {"id": 2, "value": [0.25, 0.25, 0.25]}],
        "units": "crystal",
        "name": "basis",
    },
    "lattice": {
        "type": "FCC",
        "a": 3.867,
        "b": 3.867,
        "c": 3.867,
        "alpha": 60,
        "beta": 60,
        "gamma": 60,
        "units": {"length": "angstrom", "angle": "degree"},
        "vectors": {
            "a": [3.867, 0, 0],
            "b": [1.9335000000000004, 3.348920236434424, 0],
            "c": [1.9335000000000004, 1.1163067454781415, 3.1573922784475164],
            "name": "lattice vectors",
            "alat": 1,
            "units": "angstrom",
        },
    },
    "tags": ["REST API"],
}


# ## Create material
# 
# Initialize `MaterialEndpoints` class and call `create` function to create material.

# In[ ]:


endpoint = MaterialEndpoints(*ENDPOINT_ARGS)
material = endpoint.create(CONFIG)


# ## Print new material

# In[ ]:


display_JSON(material)

