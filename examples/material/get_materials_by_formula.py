#!/usr/bin/env python
# coding: utf-8

# # Overview
# 
# Inside this example we contact [Material](https://docs.exabyte.io/api/Material/get_materials) endpoint to obtain a list materials that an account has access to. We use chemical formula to filter the list.

# # Execution
# 
# > <span style="color: orange">**NOTE**</span>: In order to run this example, an active Exabyte.io account is required. RESTful API credentials shall be updated in [settings](../settings.py). The generation of the credentials is also explained therein.
# 
# ## Import packages

# In[2]:


import json

from exabyte_api_client.endpoints.materials import MaterialEndpoints

# Import settings file
import os,sys
parent_dir = os.path.dirname(os.path.dirname(__file__))
sys.path.insert(0, parent_dir) # Insert at first entry to ensure settings.py and utils.py aren't shadowed
from settings import ENDPOINT_ARGS, ACCOUNT_ID


# ## Set Parameters
# 
# - **QUERY**: A query describing the documents to find. See [Meteor collection](https://docs.meteor.com/api/collections.html#Mongo-Collection-find) for more information.

# In[ ]:


QUERY = {
    "formula": "Si",
    "owner._id": ACCOUNT_ID
}


# ## Initialize the endpoint

# In[4]:


endpoint = MaterialEndpoints(*ENDPOINT_ARGS)


# ## List materials
# 
# Contact the endpoint to list materials according to the query above.

# In[5]:


materials = endpoint.list(QUERY)


# ## Print materials
# 
# Print the list of materials saved under the corresponding variable in pretty JSON below.

# In[6]:


print(json.dumps(materials, indent=4))

