#!/usr/bin/env python
# coding: utf-8

# <a href="https://colab.research.google.com/github/Exabyte-io/api-examples/blob/dev/examples/material/get_materials_by_formula.ipynb" target="_parent">
# <img alt="Open in Google Colab" src="https://user-images.githubusercontent.com/20477508/128780728-491fea90-9b23-495f-a091-11681150db37.jpeg" width="150" border="0">
# </a>

# # Overview
# 
# Inside this example we contact [Material](https://docs.mat3ra.com/api/Material/get_materials) endpoint to obtain a list materials that an account has access to. We use chemical formula to filter the list.

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


from utils.settings import ENDPOINT_ARGS, ACCOUNT_ID
from utils.generic import display_JSON

from exabyte_api_client.endpoints.materials import MaterialEndpoints


# ## Set Parameters
# 
# - **QUERY**: A query describing the documents to find. See [Meteor collection](https://docs.meteor.com/api/collections.html#Mongo-Collection-find) for more information.

# In[ ]:


QUERY = {"formula": "Si", "owner._id": ACCOUNT_ID}


# ## Initialize the endpoint

# In[ ]:


endpoint = MaterialEndpoints(*ENDPOINT_ARGS)


# ## List materials
# 
# Contact the endpoint to list materials according to the query above.

# In[ ]:


materials = endpoint.list(QUERY)


# ## Print materials
# 
# Print the list of materials saved under the corresponding variable in pretty JSON below.

# In[ ]:


display_JSON(materials)

