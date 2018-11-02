#!/usr/bin/env python
# coding: utf-8

# # Overview
# 
# Inside this example we contact [materials](https://docs.exabyte.io/api/#!/Materials/get_materials) endpoint to obtain a list materials that an account has access to. We use chemical formula to filter the list.

# 1. Import required packages. Adjust [settings](../settings.ipynb) as necessary.

# In[4]:


import json
import argparse
import nbimporter
nbimporter.options['only_defs'] = False

from settings import *
from endpoints.materials import MaterialEndpoints


# 2. Set QUERY in [MongoDB](https://docs.mongodb.com/manual/tutorial/query-documents/) format. 

# In[5]:


QUERY = {
    "formula": "SiGe"
}


# 3. Initialize a helper class to interact with `MaterialEndpoints`. This only has to be done once.

# In[6]:


endpoint = MaterialEndpoints(HOST, PORT, ACCOUNT_ID, AUTH_TOKEN, VERSION, SECURE)


# 4. Contact the endpoint to list materials according to the query above.

# In[7]:


materials = endpoint.list(QUERY)


# 5. Print the list of materials saved under the corresponding variable in pretty JSON below.

# In[8]:


print json.dumps(materials, indent=4)

