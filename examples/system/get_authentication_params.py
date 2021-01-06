#!/usr/bin/env python
# coding: utf-8

# # Overview
# 
# This example shows how to log in to Exabyte RESTFul API via [Login](https://docs.exabyte.io/api/API/post_login) endpoint and generate API authentication parameters.

# # Execution
# 
# > <span style="color: orange">**NOTE**</span>: In order to run this example, an active Exabyte.io account is required.
# 
# ## Import packages

# In[2]:


import json

from exabyte_api_client.endpoints.login import LoginEndpoint
from .settings import HOST, PORT, VERSION, SECURE


# ## Set Parameters
# 
# - **USERNAME**: Your Exabyte account username.
# 
# - **PASSWORD**: Your Exabyte account password.

# In[3]:


USERNAME = "YOUR_USERNANE"
PASSWORD = "YOUR_PASSWORD"


# ## Initialize the endpoint

# In[4]:


endpoint = LoginEndpoint(HOST, PORT, USERNAME, PASSWORD, VERSION, SECURE)
auth_params = endpoint.login()


# ## Print authentication parameters
# 
# Print the authentication parameters in pretty JSON below. Update [settings](../settings.py) with this parameters to be able to run other examples.

# In[5]:


print json.dumps(auth_params, indent=4)

