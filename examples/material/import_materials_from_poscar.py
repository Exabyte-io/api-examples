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

# In[]:


import os
import sys

# Import settings and utils file
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path: sys.path.append(module_path)
from settings import ENDPOINT_ARGS
from utils import ensure_packages_are_installed, display_JSON
ensure_packages_are_installed()

from exabyte_api_client.endpoints.materials import MaterialEndpoints


# ## Set Parameters
# 
# - **NAME**: material name
# - **POSCAR_PATH**: absolute path to the POSCAR file

# In[]:


NAME = "My Material"
POSCAR_PATH = "mp-978534.poscar"


# ## Import material
# 
# Initialize `MaterialEndpoints` class and call `import_from_file` function to import the material.

# In[]:


content  = ""
with open(POSCAR_PATH) as f:
    content = f.read()

endpoint = MaterialEndpoints(*ENDPOINT_ARGS)
material = endpoint.import_from_file(NAME, content)


# ## Print imported material
# 
# Print the list of imported materials in pretty JSON below.

# In[]:


display_JSON(material)

