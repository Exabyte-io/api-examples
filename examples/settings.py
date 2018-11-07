#!/usr/bin/env python
# coding: utf-8

# # Overview
# 
# This file contains settings that are used inside the examples.
# 
# Adjust settings as necessary and save the changes. When the file is saved a post-save [hook](https://jupyter-notebook.readthedocs.io/en/stable/extending/savehooks.html) is triggered to create the script (settings.py) which should be symlinked from within the example directories as notebook does not support importing from the other directories.
# 
# > <span style="color: orange">**NOTE**</span>: Jupyter may need to be reloaded on changes to load the new settings.
# 
# 
# # Settings
# 
# - **HOST**: Hostname of the RESTful API server. Defaults to platform.exabyte.io.
# 
# 
# - **PORT**: The port RESTful API server is listening on. Defaults to 443.
# 
# 
# - **VERSION**: RESTFul API version. Defaults to 2018-10-01.
# 
# 
# - **SECURE**: Whether to use secure connection. Defaults to True.
# 
# 
# - **ACCOUNT_ID**: Your account ID. See [this](../system/get_authentication_params.ipynb) example to obtain RESTful API authentication parameters.
# 
# 
# - **AUTH_TOKEN**: Your authentication token. See [this](../system/get_authentication_params.ipynb) example to obtain RESTful API authentication parameters.
# 
# 
# - **MATERIALS_PROJECT_API_KEY**: Your materials project [API](https://materialsproject.org/open) key.

# In[ ]:


PORT = 443
SECURE = True
VERSION = "2018-10-01"
HOST = "platform.exabyte.io"

ACCOUNT_ID = "YOUR_ACCOUNT_ID"
AUTH_TOKEN = "YOUR_AUTH_TOKEN"
MATERIALS_PROJECT_API_KEY = "YOUR_API_KEY"

