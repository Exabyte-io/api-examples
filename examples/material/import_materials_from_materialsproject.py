#!/usr/bin/env python
# coding: utf-8

# # Overview
# 
# This example demonstrates how to import materials from the materials project database via [Material](https://docs.exabyte.io/api/Material/post_materials_import) endpoint.

# # Execution
# 
# > <span style="color: orange">**NOTE**</span>: In order to run this example, an active Exabyte.io account is required. RESTful API credentials shall be updated in [settings](../settings.ipynb). The generation of the credentials is also explained therein.
# 
# ## Import packages

# In[2]:


import json

from endpoints.materials import MaterialEndpoints
from settings import HOST, PORT, ACCOUNT_ID, AUTH_TOKEN, VERSION, SECURE, MATERIALS_PROJECT_API_KEY


# ## Set Parameters
# 
# - **MATERIALS_PROJECT_IDS**: a list of material IDs to be imported
# 
# - **TAGS**: a list of [tags](https://docs.exabyte.io/entities-general/data/#tags) to assign to imported materials

# In[3]:


MATERIALS_PROJECT_IDS = ["mp-978534", "mp-1096549"]
TAGS = ["phase-ii", "difficulty-1"]


# The below is an embeded [iframe](https://ipython.org/ipython-doc/2/api/generated/IPython.lib.display.html) in the IPython notebook to visualize the material.

# In[4]:


from IPython.display import IFrame    
IFrame('https://materialsproject.org/materials/{}'.format(MATERIALS_PROJECT_IDS[0]), width=800, height=650)


# ## Import materials
# 
# Initialize `MaterialEndpoints` class and call `import_from_materialsproject` function to import materials.

# In[5]:


endpoint = MaterialEndpoints(HOST, PORT, ACCOUNT_ID, AUTH_TOKEN, VERSION, SECURE)
materials = endpoint.import_from_materialsproject(MATERIALS_PROJECT_API_KEY, MATERIALS_PROJECT_IDS, tags=TAGS)


# ## Print imported materials
# 
# Print the list of imported materials in pretty JSON below.

# In[6]:


print json.dumps(materials, indent=4)

