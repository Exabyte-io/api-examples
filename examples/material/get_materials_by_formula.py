#!/usr/bin/env python
# coding: utf-8

# <a href="https://colab.research.google.com/github/Exabyte-io/exabyte-api-examples/blob/feature/SOF-4618/examples/material/get_materials_by_formula.ipynb" target="_blank">Open in Google Colab</a>

# # Overview
# 
# Inside this example we contact [Material](https://docs.exabyte.io/api/Material/get_materials) endpoint to obtain a list materials that an account has access to. We use chemical formula to filter the list.

# # Complete Authorization Form and Initialize Settings
# 
# This will also determine environment and set all environment variables. We determine if we are using Jupyter Notebooks or Google Colab to run this tutorial.
# 
# ACCOUNT_ID and AUTH_TOKEN - Authentication parameters needed for when making requests to [Exabyte.io's API Endpoints](https://docs.exabyte.io/rest-api/endpoints/).
# 
# MATERIALS_PROJECT_API_KEY - Authentication parameter needed for when making requests to [Material Project's API](https://materialsproject.org/open)
# 
# ORGANIZATION_ID - Authentication parameter needed for when working with collaborative accounts https://docs.exabyte.io/collaboration/organizations/overview/
# 
# > <span style="color: orange">**NOTE**</span>: If you are running this notebook from Jupyter, the variables ACCOUNT_ID, AUTH_TOKEN, MATERIALS_PROJECT_API_KEY, and ORGANIZATION_ID should be set in the file [settings.json](../settings.json) if you need to use these variables. To obtain API token parameters, please see the following link to the documentation explaining how to get them: https://docs.exabyte.io/accounts/ui/preferences/api/



#@title Authorization Form
ACCOUNT_ID = "ACCOUNT_ID" #@param {type:"string"}
AUTH_TOKEN = "AUTH_TOKEN" #@param {type:"string"}
MATERIALS_PROJECT_API_KEY = "MATERIALS_PROJECT_API_KEY" #@param {type:"string"}
ORGANIZATION_ID  = "ORGANIZATION_ID" #@param {type:"string"}
import os, glob, sys, importlib, urllib.request

# The below execution sets up runtime using code stored remotely in a url
exec(urllib.request.urlopen('https://raw.githubusercontent.com/Exabyte-io/exabyte-api-examples/dev/examples/utils/initialize_settings.py').read())


# # Imports



from utils.generic import display_JSON
import settings; importlib.reload(settings)
from settings import ENDPOINT_ARGS, ACCOUNT_ID

from exabyte_api_client.endpoints.materials import MaterialEndpoints


# ## Set Parameters
# 
# - **QUERY**: A query describing the documents to find. See [Meteor collection](https://docs.meteor.com/api/collections.html#Mongo-Collection-find) for more information.



QUERY = {
    "formula": "Si",
    "owner._id": ACCOUNT_ID
}


# ## Initialize the endpoint



endpoint = MaterialEndpoints(*ENDPOINT_ARGS)


# ## List materials
# 
# Contact the endpoint to list materials according to the query above.



materials = endpoint.list(QUERY)


# ## Print materials
# 
# Print the list of materials saved under the corresponding variable in pretty JSON below.



display_JSON(materials)
