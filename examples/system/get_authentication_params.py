#!/usr/bin/env python
# coding: utf-8

# # Overview
# 
# This example shows how to log in to Exabyte RESTFul API via [Login](https://docs.exabyte.io/api/API/post_login) endpoint and generate API authentication parameters.

# 1. Import required packages. Adjust [settings](../settings.ipynb) as necessary.

# In[1]:


import json
import argparse

from settings import *
from endpoints.login import LoginEndpoint


# 2. Set Parameters
# 
#     - **USERNAME**: Your Exabyte account username.
# 
#     - **PASSWORD**: Your Exabyte account password.

# In[2]:


USERNAME = "YOUR_USERNANE"
PASSWORD = "YOUR_PASSWORD"


# 3. Initialize `LoginEndpoint` and call `login` function to generate authentication parameters.

# In[ ]:


endpoint = LoginEndpoint(HOST, PORT, USERNAME, PASSWORD, VERSION, SECURE)
auth_params = endpoint.login()


# 4. Print the authentication parameters in pretty JSON below. Save them in a secure place as they are required in other calls to the API.

# In[4]:


print json.dumps(auth_params, indent=4)

