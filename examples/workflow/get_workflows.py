#!/usr/bin/env python
# coding: utf-8

# # Overview
# 
# Inside this example we contact [Workflow](https://docs.exabyte.io/api/Workflows/get_workflows) endpoint to obtain a list of workflows that an account has access to.

# # Execution
# 
# > <span style="color: orange">**NOTE**</span>: In order to run this example, an active Exabyte.io account is required. RESTful API credentials shall be updated in [settings](../settings.py). The generation of the credentials is also explained therein.
# 
# ## Import packages

# In[]:


import os
import sys

# Import settings and utils file
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path: sys.path.append(module_path)
from settings import ENDPOINT_ARGS, ACCOUNT_ID
from utils import ensure_packages_are_installed, display_JSON
ensure_packages_are_installed()

from exabyte_api_client.endpoints.workflows import WorkflowEndpoints


# ## Set Parameters
# 
# - **QUERY**: A query describing the documents to find. See [Meteor collection](https://docs.meteor.com/api/collections.html#Mongo-Collection-find) for more information. 
# 
# - **limit**: Maximum number of results to return. See [Meteor collection](https://docs.meteor.com/api/collections.html#Mongo-Collection-find) for more information.

# In[]:


QUERY = {
    "name": "Total Energy",
    "owner._id": ACCOUNT_ID
}

OPTIONS = {
    "limit": 2
}


# ## Initialize the endpoint
# 
# Initialize a helper class to interact with `WorkflowEndpoints`. This only has to be done once.

# In[]:


endpoint = WorkflowEndpoints(*ENDPOINT_ARGS)


# ## List workflows
# 
# Contact the endpoint to list workflows according to the query above.

# In[]:


workflows = endpoint.list(QUERY, OPTIONS)


# ## Print workflows
# 
# Print the list of workflows saved under the corresponding variable in pretty JSON below.

# In[]:


display_JSON(workflows)

