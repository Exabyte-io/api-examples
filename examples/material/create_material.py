#!/usr/bin/env python
# coding: utf-8

# [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Exabyte-io/exabyte-api-examples/blob/feature/SOF-4400-skinny-req/examples/material/create_material.ipynb)

# # Overview
# 
# In this example we create a material from a JSON config with [tags](https://docs.exabyte.io/entities-general/data/#tags) to identify the material.

# # Complete Authorization Form and Initialize Settings
# 
# This will also determine environment and set all environment variables. We determine if we are using Jupyter Notebooks or Google Colab to run this tutorial.



#@title Authorization Form
ACCOUNT_ID = "ACCOUNT_ID" #@param {type:"string"}
AUTH_TOKEN = "AUTH_TOKEN" #@param {type:"string"}
MATERIALS_PROJECT_API_KEY = "MATERIALS_PROJECT_API_KEY" #@param {type:"string"}
ORGANIZATION_ID  = "ORGANIZATION_ID" #@param {type:"string"}
import os, glob, sys, importlib, urllib.request

# The below execution sets up runtime using code stored remotely in a url
exec(urllib.request.urlopen('https://raw.githubusercontent.com/Exabyte-io/exabyte-api-examples/feature/SOF-4400-skinny-req/examples/utils/notebooks-executables.py').read())


# # Imports



from utils.generic import display_JSON
import settings; importlib.reload(settings); from settings import ENDPOINT_ARGS
from exabyte_api_client.endpoints.materials import MaterialEndpoints


# ## Create material config
# 
# Create material config in JSON format. See [Material](https://docs.exabyte.io/api/Material/put_materials_create) endpoint for more information about material config format.



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



endpoint = MaterialEndpoints(*ENDPOINT_ARGS)
material = endpoint.create(CONFIG)


# ## Print new material



display_JSON(material)

