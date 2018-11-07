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

# In[4]:


import json

from endpoints.materials import MaterialEndpoints
from settings import HOST, PORT, ACCOUNT_ID, AUTH_TOKEN, VERSION, SECURE


# ## Set QUERY in [MongoDB](https://docs.mongodb.com/manual/tutorial/query-documents/) format

# In[5]:


QUERY = {
    "formula": "SiGe",
    "owner.slug": "demo"
}


# ## Initialize the endpoint

# In[6]:


endpoint = MaterialEndpoints(HOST, PORT, ACCOUNT_ID, AUTH_TOKEN, VERSION, SECURE)


# ## List materials
# 
# Contact the endpoint to list materials according to the query above.

# In[7]:


materials = endpoint.list(QUERY)


# ## Print materials
# 
# Print the list of materials saved under the corresponding variable in pretty JSON below.

# In[8]:


print json.dumps(materials, indent=4)

