#!/usr/bin/env python
# coding: utf-8

# [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Exabyte-io/exabyte-api-examples/blob/feature/SOF-4400-skinny-req/examples/material/create_material.ipynb)

# # Overview
# 
# In this example we create a material from a JSON config with [tags](https://docs.exabyte.io/entities-general/data/#tags) to identify the material.

# # Complete Authorization Form and Initialize Settings
# 
# This will also determine environment and set all environment variables. We determine if we are using Jupyter Notebooks or Google Colab to run this tutorial. If google drive, we will have some extra steps to do.

# In[ ]:


#@title Authorization Form
ACCOUNT_ID = "" #@param {type:"string"}
AUTH_TOKEN = "" #@param {type:"string"}
MATERIALS_PROJECT_API_KEY = "" #@param {type:"string"}

import os

if 'google.colab' in str(get_ipython()):
  print('Running on CoLab')
  os.environ['notebook_environment'] = "Colab"
elif 'ZMQInteractiveShell'  in str(get_ipython()):
  print('Running in Jupyter')
  os.environ['notebook_environment'] = "Jupyter"
else:
  print('Unknown Environment')
  os.environ['notebook_environment'] = ""
    
if os.environ['notebook_environment'] == "Colab":
    get_ipython().system('git clone -b feature/SOF-4400-skinny-req https://github.com/Exabyte-io/exabyte-api-examples.git')
    get_ipython().run_line_magic('cd', '/content/exabyte-api-examples/examples/material')
    get_ipython().system('pip install --no-deps -r ../../requirements-colab.txt')

with open('../settings.py') as settings:
    settings_filelines = settings.readlines()
with open('../settings.py', "w") as settings:
    for line in settings_filelines:
        if 'ACCOUNT_ID = ' in line:
            settings.write('%s "%s"\n' % ("ACCOUNT_ID =", ACCOUNT_ID))
        elif 'AUTH_TOKEN = ' in line:
            settings.write('%s "%s"\n' % ("AUTH_TOKEN =", AUTH_TOKEN))
        elif 'MATERIALS_PROJECT_API_KEY = ' in line:
            settings.write('%s "%s"\n' % ("MATERIALS_PROJECT_API_KEY =", MATERIALS_PROJECT_API_KEY))
        else:
            settings.write(line)


# # Imports

# In[ ]:


import sys
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path: sys.path.append(module_path)
from settings import ENDPOINT_ARGS
from utils.generic import display_JSON

if os.environ['notebook_environment'] == "Jupyter":
    from utils.generic import ensure_packages_are_installed
    ensure_packages_are_installed()

from exabyte_api_client.endpoints.materials import MaterialEndpoints


# ## Create material config
# 
# Create material config in JSON format. See [Material](https://docs.exabyte.io/api/Material/put_materials_create) endpoint for more information about material config format.

# In[ ]:


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

# In[ ]:


endpoint = MaterialEndpoints(*ENDPOINT_ARGS)
material = endpoint.create(CONFIG)


# ## Print new material

# In[ ]:


display_JSON(material)


# In[ ]:




