#!/usr/bin/env python
# coding: utf-8

# <a href="https://colab.research.google.com/github/Exabyte-io/exabyte-api-examples/blob/feature/SOF-4685/examples/material/import_materials_from_materialsproject.ipynb" target="_blank">Open in Google Colab</a>

# # Overview
# 
# This example demonstrates how to import materials from the materials project database via [Material](https://docs.exabyte.io/api/Material/post_materials_import) endpoint.

# # Complete Authorization Form and Initialize Settings
# 
# This will also determine environment and set all environment variables. We determine if we are using Jupyter Notebooks or Google Colab to run this tutorial.
# 
# If you are running this notebook from Google Colab, Colab takes ~1 min to execute the following cell.
# 
# ACCOUNT_ID and AUTH_TOKEN - Authentication parameters needed for when making requests to [Exabyte.io's API Endpoints](https://docs.exabyte.io/rest-api/endpoints/).
# 
# MATERIALS_PROJECT_API_KEY - Authentication parameter needed for when making requests to [Material Project's API](https://materialsproject.org/open)
# 
# ORGANIZATION_ID - Authentication parameter needed for when working with collaborative accounts https://docs.exabyte.io/collaboration/organizations/overview/
# 
# > <span style="color: orange">**NOTE**</span>: If you are running this notebook from Jupyter, the variables ACCOUNT_ID, AUTH_TOKEN, MATERIALS_PROJECT_API_KEY, and ORGANIZATION_ID should be set in the file [settings.json](../settings.json) if you need to use these variables. To obtain API token parameters, please see the following link to the documentation explaining how to get them: https://docs.exabyte.io/accounts/ui/preferences/api/

# In[]:


#@title Authorization Form
ACCOUNT_ID = "ACCOUNT_ID" #@param {type:"string"}
AUTH_TOKEN = "AUTH_TOKEN" #@param {type:"string"}
MATERIALS_PROJECT_API_KEY = "MATERIALS_PROJECT_API_KEY" #@param {type:"string"}
ORGANIZATION_ID  = "ORGANIZATION_ID" #@param {type:"string"}
import os, glob, sys, importlib, urllib.request

# The below execution sets up runtime using code stored remotely in a url
if 'isExecuted' not in globals():
    exec(urllib.request.urlopen('https://raw.githubusercontent.com/Exabyte-io/exabyte-api-examples/dev/examples/utils/initialize_settings.py').read())
    isExecuted = True


# # Imports

# In[]:


from utils.generic import display_JSON
import settings; importlib.reload(settings)
from settings import ENDPOINT_ARGS, MATERIALS_PROJECT_API_KEY

from exabyte_api_client.endpoints.materials import MaterialEndpoints


# ## Set Parameters
# 
# - **MATERIALS_PROJECT_IDS**: a list of material IDs to be imported
# 
# - **TAGS**: a list of [tags](https://docs.exabyte.io/entities-general/data/#tags) to assign to imported materials

# In[]:


MATERIALS_PROJECT_IDS = ["mp-978534", "mp-1096549"]
TAGS = ["tag1", "tag2"]


# ## Import materials
# 
# Initialize `MaterialEndpoints` class and call `import_from_materialsproject` function to import materials.

# In[]:


endpoint = MaterialEndpoints(*ENDPOINT_ARGS)
materials = endpoint.import_from_materialsproject(MATERIALS_PROJECT_API_KEY, MATERIALS_PROJECT_IDS, tags=TAGS)


# ## Print imported materials
# 
# Print the list of imported materials in pretty JSON below.

# In[]:


display_JSON(materials)
