#!/usr/bin/env python
# coding: utf-8

# # Settings
# 
# This file contains common variables that are used inside the examples.
# 
# Adjust settings as necessary and save the changes. When the file is saved a post-save [hook](https://jupyter-notebook.readthedocs.io/en/stable/extending/savehooks.html) is triggered to create the script (settings.py). The latter is linked into each of the example directories to facilitate imports.
# 
# > <span style="color: orange">**NOTE**:</span> Jupyter may need to be reloaded on changes to load the new settings.
# 
# 
# ## Variables
# 
# ### Account settings. 
# 
# Allow for the exeucution of the API calls under a particular account. Need a one-time adjustments for examples to work.
# 
# - **ACCOUNT_SLUG**: Computer-friendly name (or "slug") of the account (eg. "exabyte-io). See [main documentation](https://docs.exabyte.io) for more explanation.
# - **ACCOUNT_ID**: Account ID string. See [this](../system/get_authentication_params.ipynb) example to obtain RESTful API authentication parameters.
# - **AUTH_TOKEN**: Account authentication token. See [this](../system/get_authentication_params.ipynb) example to obtain RESTful API authentication parameters.
# 
# - **MATERIALS_PROJECT_API_KEY**: Your materials project [API](https://materialsproject.org/open) key.
# 
# ### Advanced settings. 
# 
# Facilitate connectivity with the API. Should not need adjustments.
# 
# - **HOST**: Hostname of the RESTful API server. Defaults to platform.exabyte.io.
# - **PORT**: The port RESTful API server is listening on. Defaults to 443.
# - **VERSION**: RESTFul API version. Defaults to 2018-10-01.
# - **SECURE**: Whether to use secure connection. Defaults to True.

# In[ ]:


# Account settings. Need a one-time adjustments for examples to work.

ACCOUNT_SLUG = "ACCOUNT_SLUG"

ACCOUNT_ID = "ACCOUNT_ID"
AUTH_TOKEN = "AUTH_TOKEN"

MATERIALS_PROJECT_API_KEY = "MATERIALS_PROJECT_API_KEY"

# Advanced settings. Should not need adjustments

PORT = 443
SECURE = True
VERSION = "2018-10-01"
HOST = "platform.exabyte.io"
ENDPOINT_ARGS = [HOST, PORT, ACCOUNT_ID, AUTH_TOKEN, VERSION, SECURE]

