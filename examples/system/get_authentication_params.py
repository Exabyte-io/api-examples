#!/usr/bin/env python
# coding: utf-8

# <a href="https://colab.research.google.com/github/Exabyte-io/exabyte-api-examples/blob/feature/SOF-4685/examples/system/get_authentication_params.ipynb" target="_parent">
# <img alt="Open in Google Colab" src="https://user-images.githubusercontent.com/20477508/128780728-491fea90-9b23-495f-a091-11681150db37.jpeg" width="150" border="0">
# </a>

# # Overview
# 
# This example shows how to log in to Exabyte RESTFul API via [Login](https://docs.exabyte.io/api/API/post_login) endpoint and generate API authentication parameters.
# 
# Here, we execute a remote URL to set our notebook environment. Do not edit the following cell's contents.

# In[]:


ACCOUNT_ID = AUTH_TOKEN = MATERIALS_PROJECT_API_KEY = ORGANIZATION_ID = ''
import os, glob, sys, importlib, urllib.request

# The below execution sets up runtime using code stored remotely in a url
exec(urllib.request.urlopen('https://raw.githubusercontent.com/Exabyte-io/exabyte-api-examples/feature/SOF-4685/examples/utils/initialize_settings.py').read())


# # Execution
# 
# > <span style="color: orange">**NOTE**</span>: In order to run this example, an active Exabyte.io account is required.

# ## Import packages

# In[]:


from settings import HOST, PORT, VERSION, SECURE
from utils.generic import display_JSON

from exabyte_api_client.endpoints.login import LoginEndpoint


# ## Set Parameters
# 
# - **USERNAME**: Your Exabyte account username.
# 
# - **PASSWORD**: Your Exabyte account password.

# In[]:


USERNAME = "YOUR_USERNANE"
PASSWORD = "YOUR_PASSWORD"


# ## Initialize the endpoint

# In[]:


endpoint = LoginEndpoint(HOST, PORT, USERNAME, PASSWORD, VERSION, SECURE)
auth_params = endpoint.login()


# ## Print authentication parameters
# 
# Print the authentication parameters in pretty JSON below. Update [settings](../settings.py) with this parameters to be able to run other examples.

# In[]:


display_JSON(auth_params)
