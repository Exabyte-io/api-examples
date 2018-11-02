#!/usr/bin/env python
# coding: utf-8

# # Overview
# 
# This example demonstrates how to import materials from the materials project database via [Material](https://docs.exabyte.io/api/Material/post_materials_import) endpoint.

# 1. Import required packages. Adjust [settings](../settings.ipynb) as necessary.

# In[2]:


import json
import argparse
import nbimporter
nbimporter.options['only_defs'] = False

from settings import *
from endpoints.materials import MaterialEndpoints


# 2. Set Parameters
# 
#     - **MATERIALS_PROJECT_API_KEY**: Your materials project API key.
#     
#     - **MATERIALS_PROJECT_IDS**: A list of material IDs you would like to import.
#     
#     - **TAGS**: A list of tags you want to assign to imported materials.

# In[2]:


MATERIALS_PROJECT_API_KEY = "YOUR_API_KEY"
MATERIALS_PROJECT_IDS = ["mp-978534", "mp-1096549"]
TAGS = ["phase-ii", "difficulty-1"]


# The below is an embeded [iframe](https://ipython.org/ipython-doc/2/api/generated/IPython.lib.display.html) in the IPython notebook to visualize the material.

# In[8]:


from IPython.display import IFrame    
IFrame('https://materialsproject.org/materials/{}'.format(MATERIALS_PROJECT_IDS[0]), width=800, height=650)


# 3. Initialize `MaterialEndpoints` class and call `import_from_materialsproject` function to import materials.

# In[7]:


endpoint = MaterialEndpoints(HOST, PORT, ACCOUNT_ID, AUTH_TOKEN, VERSION, SECURE)
materials = endpoint.import_from_materialsproject(MATERIALS_PROJECT_API_KEY, MATERIALS_PROJECT_IDS, tags=TAGS)


# 4. Print the list of imported materials in pretty JSON below.

# In[6]:


print json.dumps(materials, indent=4)

