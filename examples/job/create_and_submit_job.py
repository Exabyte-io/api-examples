#!/usr/bin/env python
# coding: utf-8

# <a href="https://colab.research.google.com/github/Exabyte-io/exabyte-api-examples/blob/feature/SOF-4618/examples/job/create_and_submit_job.ipynb" target="_blank">Open in Google Colab</a>

# # Overview
# 
# This example demonstrates how to create and submit a job via [Job](https://docs.exabyte.io/api/Job/put_jobs_create) endpoints.

# # Complete Authorization Form and Initialize Settings
# 
# This will also determine environment and set all environment variables. We determine if we are using Jupyter Notebooks or Google Colab to run this tutorial.
# 
# ACCOUNT_ID and AUTH_TOKEN - Authentication parameters needed for when making requests to [Exabyte.io's API Endpoints](https://docs.exabyte.io/rest-api/endpoints/).
# 
# MATERIALS_PROJECT_API_KEY - Authentication parameter needed for when making requests to [Material Project's API](https://materialsproject.org/open)
# 
# ORGANIZATION_ID - Authentication parameter needed for when working with collaborative accounts https://docs.exabyte.io/collaboration/organizations/overview/
# 
# > <span style="color: orange">**NOTE**</span>: If you are running this notebook from Jupyter, the variables ACCOUNT_ID, AUTH_TOKEN, MATERIALS_PROJECT_API_KEY, and ORGANIZATION_ID should be set in the file [settings.json](../settings.json) if you need to use these variables. To obtain API token parameters, please see the following link to the documentation explaining how to get them: https://docs.exabyte.io/accounts/ui/preferences/api/



#@title Authorization Form
ACCOUNT_ID = "ACCOUNT_ID" #@param {type:"string"}
AUTH_TOKEN = "AUTH_TOKEN" #@param {type:"string"}
MATERIALS_PROJECT_API_KEY = "MATERIALS_PROJECT_API_KEY" #@param {type:"string"}
ORGANIZATION_ID  = "ORGANIZATION_ID" #@param {type:"string"}
import os, glob, sys, importlib, urllib.request

# The below execution sets up runtime using code stored remotely in a url
exec(urllib.request.urlopen('https://raw.githubusercontent.com/Exabyte-io/exabyte-api-examples/dev/examples/utils/initialize_settings.py').read())




import settings; importlib.reload(settings)
from settings import ENDPOINT_ARGS, ACCOUNT_ID
from utils.generic import display_JSON

from exabyte_api_client.endpoints.jobs import JobEndpoints
from exabyte_api_client.endpoints.materials import MaterialEndpoints
from exabyte_api_client.endpoints.workflows import WorkflowEndpoints


# ## Initialize the endpoints



job_endpoints = JobEndpoints(*ENDPOINT_ARGS)
material_endpoints = MaterialEndpoints(*ENDPOINT_ARGS)
workflow_endpoints = WorkflowEndpoints(*ENDPOINT_ARGS)


# Set job name.



JOB_NAME = "TEST JOB"


# ## Retrieve IDs
# 
# Default account's materail and workflow are used in this example to create the job. Adjust the queries to use different material and workflow.



default_material = material_endpoints.list({"isDefault": True, "owner._id": ACCOUNT_ID})[0]
default_workflow = workflow_endpoints.list({"isDefault": True, "owner._id": ACCOUNT_ID})[0]

material_id = default_material["_id"]
workflow_id = default_workflow["_id"]
owner_id = default_material["owner"]["_id"]


# ## Create job config
# 
# The job belongs to user's default account and it is created inside the defauult account's project. 



config = {
    "owner": {
        "_id": owner_id
    },
    "_material": {
        "_id": material_id
    },
    "workflow": {
        "_id": workflow_id
    },
    "name": JOB_NAME
}


# ## Create and submit job



job = job_endpoints.create(config)
job_endpoints.submit(job['_id'])


# ## Print the job
# 
# Print the job in pretty JSON below. Check `status` field to make sure job is submiited.



job = job_endpoints.get(job['_id'])
display_JSON(job)
