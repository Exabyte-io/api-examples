#!/usr/bin/env python
# coding: utf-8

# # Overview
# 
# This example demonstrates how to import materials from the materials project database via [Material](https://docs.exabyte.io/api/Material/post_materials_import) endpoint.

# # Execution
# 
# > <span style="color: orange">**NOTE**</span>: In order to run this example, an active Exabyte.io account is required. RESTful API credentials shall be updated in [settings](../settings.py). The generation of the credentials is also explained therein.
# 
# ## Import packages

# In[9]:


import os
import sys
from IPython.display import JSON

# Import settings file and utils
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path: sys.path.append(module_path)
from settings import ENDPOINT_ARGS, MATERIALS_PROJECT_API_KEY
from utils import ensure_installed

ensure_installed("exabyte_api_client")
from exabyte_api_client.endpoints.materials import MaterialEndpoints


# ## Set Parameters
# 
# - **MATERIALS_PROJECT_IDS**: a list of material IDs to be imported
# 
# - **TAGS**: a list of [tags](https://docs.exabyte.io/entities-general/data/#tags) to assign to imported materials

# In[10]:


MATERIALS_PROJECT_IDS = ["mp-978534", "mp-1096549"]
TAGS = ["tag1", "tag2"]


# ## Import materials
# 
# Initialize `MaterialEndpoints` class and call `import_from_materialsproject` function to import materials.

# In[11]:


endpoint = MaterialEndpoints(*ENDPOINT_ARGS)
materials = endpoint.import_from_materialsproject(MATERIALS_PROJECT_API_KEY, MATERIALS_PROJECT_IDS, tags=TAGS)


# ## Print imported materials
# 
# Print the list of imported materials in pretty JSON below.

# In[12]:


JSON(materials)

