#!/usr/bin/env python
# coding: utf-8

# # Overview
# 
# In this example we create a material from a JSON config with certain tags to identify the material.

# 1. Import required packages. Adjust [settings](../settings.ipynb) as necessary.

# In[2]:


import json
import argparse

from settings import *
from endpoints.materials import MaterialEndpoints


# 2. Create material config in JSON format. See [Material](https://docs.exabyte.io/api/Material/put_materials_create) endpoint for more information about material config format.

# In[3]:


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


# 3. Initialize `MaterialEndpoints` class and call `create` function to create material.

# In[4]:


endpoint = MaterialEndpoints(HOST, PORT, ACCOUNT_ID, AUTH_TOKEN, VERSION, SECURE)
material = endpoint.create(CONFIG)


# 4. Print new material in pretty JSON below.

# In[5]:


print json.dumps(material, indent=4)

