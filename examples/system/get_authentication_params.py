#!/usr/bin/env python
# coding: utf-8

# <a href="https://colab.research.google.com/github/Exabyte-io/api-examples/blob/dev/examples/system/get_authentication_params.ipynb" target="_parent">
# <img alt="Open in Google Colab" src="https://user-images.githubusercontent.com/20477508/128780728-491fea90-9b23-495f-a091-11681150db37.jpeg" width="150" border="0">
# </a>

# # Overview
# 
# This example shows how to log in to Mat3ra RESTFul API via [Login](https://docs.mat3ra.com/api/API/post_login) endpoint and generate API authentication parameters. Note that it is also possible to generate an API token manually via platform web UI as described in the project [README](https://github.com/Exabyte-io/api-examples/tree/dev#usage).

# # Execution
# 
# > <span style="color: orange">**NOTE**</span>: In order to run this example, an active Mat3ra.com account is required.

# ## Set Parameters
# 
# - **USERNAME**: Your Mat3ra account username.
# 
# - **PASSWORD**: Your Mat3ra account password.

# In[ ]:


# @title Authorization Form
USERNAME = "YOUR_USERNAME"  # @param {type:"string"}

# avoid storing password in plaintext
from getpass import getpass

PASSWORD = getpass("Please enter password: ")

import os

if "COLAB_JUPYTER_IP" in os.environ:
    get_ipython().system('GIT_BRANCH="dev"; export GIT_BRANCH; curl -s "https://raw.githubusercontent.com/Exabyte-io/api-examples/${GIT_BRANCH}/scripts/env.sh" | bash')


# ## Import packages

# In[ ]:


from utils.settings import HOST, PORT, VERSION, SECURE
from utils.generic import display_JSON

from exabyte_api_client.endpoints.login import LoginEndpoint


# ## Initialize the endpoint

# In[ ]:


endpoint = LoginEndpoint(HOST, PORT, USERNAME, PASSWORD, VERSION, SECURE)
auth_params = endpoint.login()


# ## Print authentication parameters
# 
# Print the authentication parameters in pretty JSON below. Update [settings](../../utils/settings.json) with this parameters to be able to run other examples.

# In[ ]:


display_JSON(auth_params)

