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

# In[1]:


import os
import sys
from IPython.display import JSON

# Install Pandas if it isn't present
try:
    import pandas as pd
except ModuleNotFoundError:
    import subprocess, sys
    subprocess.call([sys.executable,'-m','pip','install','pandas==1.1.4'])
    import pandas as pd
# Install Tabulate if it isn't present
try:
    import tabulate
except ModuleNotFoundError:
    import subprocess, sys
    subprocess.call([sys.executable,'-m','pip','install','tabulate==0.8.2'])
    import tabulate
# Install the API Client if it isn't present
try:
    import exabyte_api_client
except ModuleNotFoundError:
    import subprocess, sys
    subprocess.call([sys.executable,'-m','pip','install','exabyte_api_client==2020.10.19'])
from exabyte_api_client.endpoints.materials import MaterialEndpoints

# Import settings file
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path: sys.path.append(module_path)
from settings import ENDPOINT_ARGS, ACCOUNT_ID


# ## Set Parameters
# 
# - **QUERY**: A query describing the documents to find. See [Meteor collection](https://docs.meteor.com/api/collections.html#Mongo-Collection-find) for more information.

# In[2]:


QUERY = {
    "formula": "Si",
    "owner._id": ACCOUNT_ID
}


# ## Initialize the endpoint

# In[3]:


endpoint = MaterialEndpoints(*ENDPOINT_ARGS)


# ## List materials
# 
# Contact the endpoint to list materials according to the query above.

# In[4]:


materials = endpoint.list(QUERY)


# ## Print materials
# 
# Print the list of materials saved under the corresponding variable in pretty JSON below.

# In[8]:


JSON(materials)


# In[ ]:




