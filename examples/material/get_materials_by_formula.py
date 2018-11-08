#!/usr/bin/env python
# coding: utf-8

# # Overview
# 
# Inside this example we contact [Material](https://docs.exabyte.io/api/Material/get_materials) endpoint to obtain a list materials that an account has access to. We use chemical formula to filter the list.

# # Execution
# 
# > <span style="color: orange">**NOTE**</span>: In order to run this example, an active Exabyte.io account is required. RESTful API credentials shall be updated in [settings](../settings.ipynb). The generation of the credentials is also explained therein.
# 
# ## Import packages

# In[2]:


import json

from settings import ENDPOINT_ARGS
from endpoints.materials import MaterialEndpoints


# ## Set Parameters
# 
# - **QUERY**: A query describing the documents to find. See [Meteor collection](https://docs.meteor.com/api/collections.html#Mongo-Collection-find) for more information.

# In[3]:


QUERY = {
    "formula": "Si",
    "owner.slug": "demo"
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


print json.dumps(materials, indent=4)

