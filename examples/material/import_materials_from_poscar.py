#!/usr/bin/env python
# coding: utf-8

# <a href="https://colab.research.google.com/github/Exabyte-io/exabyte-api-examples/blob/feature/SOF-4618/examples/material/import_materials_from_poscar.ipynb" target="_blank">Open in Google Colab</a>

# # Overview
# 
# This example demonstrates how to import a material from a POSCAR file via [Material](https://docs.exabyte.io/api/Material/post_materials_import) endpoints.

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
from settings import ENDPOINT_ARGS

from exabyte_api_client.endpoints.materials import MaterialEndpoints


# ## Set Parameters
# 
# - **NAME**: material name
# - **POSCAR_PATH**: absolute path to the POSCAR file



NAME = "My Material"
POSCAR_PATH = "../assets/mp-978534.poscar"


# ## Import material
# 
# Initialize `MaterialEndpoints` class and call `import_from_file` function to import the material.



content  = ""
with open(POSCAR_PATH) as f:
    content = f.read()

endpoint = MaterialEndpoints(*ENDPOINT_ARGS)
material = endpoint.import_from_file(NAME, content)


# ## Print imported material
# 
# Print the list of imported materials in pretty JSON below.



display_JSON(material)
