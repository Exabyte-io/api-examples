#!/usr/bin/env python
# coding: utf-8

# # Overview
# 
# In this example we create a material from a JSON config with [tags](https://docs.exabyte.io/entities-general/data/#tags) to identify the material.

# # Execution
# 
# > <span style="color: orange">**NOTE**</span>: In order to run this example, an active Exabyte.io account is required. RESTful API credentials shall be updated in [settings](../settings.py). The generation of the credentials is also explained therein.
# 
# ## Import packages

# In[1]:


import json
import os
import sys

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
from settings import ENDPOINT_ARGS


# ## Create material config
# 
# Create material config in JSON format. See [Material](https://docs.exabyte.io/api/Material/put_materials_create) endpoint for more information about material config format.

# In[2]:


CONFIG = {
    "name": "TEST MATERIAL",
    "basis": {
        "elements": [
            {
                "id": 1,
                "value": "Si"
            },
            {
                "id": 2,
                "value": "Si"
            }
        ],
        "coordinates": [
            {
                "id": 1,
                "value": [
                    0,
                    0,
                    0
                ]
            },
            {
                "id": 2,
                "value": [
                    0.25,
                    0.25,
                    0.25
                ]
            }
        ],
        "units": "crystal",
        "name": "basis"
    },
    "lattice": {
        "type": "FCC",
        "a": 3.867,
        "b": 3.867,
        "c": 3.867,
        "alpha": 60,
        "beta": 60,
        "gamma": 60,
        "units": {
            "length": "angstrom",
            "angle": "degree"
        },
        "vectors": {
            "a": [
                3.867,
                0,
                0
            ],
            "b": [
                1.9335000000000004,
                3.348920236434424,
                0
            ],
            "c": [
                1.9335000000000004,
                1.1163067454781415,
                3.1573922784475164
            ],
            "name": "lattice vectors",
            "alat": 1,
            "units": "angstrom"
        }
    },
    "tags": [
        "REST API"
    ]
}


# ## Create material
# 
# Initialize `MaterialEndpoints` class and call `create` function to create material.

# In[3]:


endpoint = MaterialEndpoints(*ENDPOINT_ARGS)
material = endpoint.create(CONFIG)


# ## Print new material

# In[4]:


print(json.dumps(material, indent=4))

