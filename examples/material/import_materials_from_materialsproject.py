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

# In[2]:


import json

from exabyte_api_client.endpoints.materials import MaterialEndpoints
from .settings import ENDPOINT_ARGS, MATERIALS_PROJECT_API_KEY


# ## Set Parameters
# 
# - **MATERIALS_PROJECT_IDS**: a list of material IDs to be imported
# 
# - **TAGS**: a list of [tags](https://docs.exabyte.io/entities-general/data/#tags) to assign to imported materials

# In[3]:


MATERIALS_PROJECT_IDS = ["mp-978534", "mp-1096549"]
TAGS = ["tag1", "tag2"]


# ## Import materials
# 
# Initialize `MaterialEndpoints` class and call `import_from_materialsproject` function to import materials.

# In[5]:


endpoint = MaterialEndpoints(*ENDPOINT_ARGS)
materials = endpoint.import_from_materialsproject(MATERIALS_PROJECT_API_KEY, MATERIALS_PROJECT_IDS, tags=TAGS)


# ## Print imported materials
# 
# Print the list of imported materials in pretty JSON below.

# In[6]:


print(json.dumps(materials, indent=4))

