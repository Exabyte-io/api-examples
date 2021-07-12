#!/usr/bin/env python
# coding: utf-8

# <a href="https://colab.research.google.com/github/Exabyte-io/exabyte-api-examples/blob/feature/SOF-4618/examples/system/get_authentication_params.ipynb" target="_blank">Open in Google Colab</a>

# # Overview
# 
# This example shows how to log in to Exabyte RESTFul API via [Login](https://docs.exabyte.io/api/API/post_login) endpoint and generate API authentication parameters.
# 
# Here, we execute a remote URL to set our notebook environment. Do not edit the following cell's contents.



ACCOUNT_ID = AUTH_TOKEN = MATERIALS_PROJECT_API_KEY = ORGANIZATION_ID = ''
import os, glob, sys, importlib, urllib.request

# The below execution sets up runtime using code stored remotely in a url
exec(urllib.request.urlopen('https://raw.githubusercontent.com/Exabyte-io/exabyte-api-examples/dev/examples/utils/initialize_settings.py').read())


# # Execution
# 
# > <span style="color: orange">**NOTE**</span>: In order to run this example, an active Exabyte.io account is required.

# ## Import packages



from settings import HOST, PORT, VERSION, SECURE
from utils.generic import display_JSON

from exabyte_api_client.endpoints.login import LoginEndpoint


# ## Set Parameters
# 
# - **USERNAME**: Your Exabyte account username.
# 
# - **PASSWORD**: Your Exabyte account password.



USERNAME = "YOUR_USERNANE"
PASSWORD = "YOUR_PASSWORD"


# ## Initialize the endpoint



endpoint = LoginEndpoint(HOST, PORT, USERNAME, PASSWORD, VERSION, SECURE)
auth_params = endpoint.login()


# ## Print authentication parameters
# 
# Print the authentication parameters in pretty JSON below. Update [settings](../settings.py) with this parameters to be able to run other examples.



display_JSON(auth_params)
