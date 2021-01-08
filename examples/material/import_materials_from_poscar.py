#!/usr/bin/env python
# coding: utf-8

# # Overview
# 
# This example demonstrates how to import a material from a POSCAR file via [Material](https://docs.exabyte.io/api/Material/post_materials_import) endpoints.

# # Execution
# 
# > <span style="color: orange">**NOTE**</span>: In order to run this example, an active Exabyte.io account is required. RESTful API credentials shall be updated in [settings](../settings.py). The generation of the credentials is also explained therein.
# 
# ## Import packages

# In[4]:


import json

from exabyte_api_client.endpoints.materials import MaterialEndpoints

# Import settings file
import os,sys
parent_dir = os.path.dirname(os.path.dirname(__file__))
sys.path.insert(0, parent_dir) # Insert at first entry to ensure settings.py and utils.py aren't shadowed
from settings import ENDPOINT_ARGS


# ## Set Parameters
# 
# - **NAME**: material name
# - **POSCAR_PATH**: absolute path to the POSCAR file

# In[5]:


NAME = "My Material"
POSCAR_PATH = "mp-978534.poscar"


# ## Import material
# 
# Initialize `MaterialEndpoints` class and call `import_from_file` function to import the material.

# In[6]:


content  = ""
with open(POSCAR_PATH) as f:
    content = f.read()

endpoint = MaterialEndpoints(*ENDPOINT_ARGS)
material = endpoint.import_from_file(NAME, content)


# ## Print imported material
# 
# Print the list of imported materials in pretty JSON below.

# In[7]:


print(json.dumps(material, indent=4))

