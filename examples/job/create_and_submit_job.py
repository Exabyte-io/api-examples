#!/usr/bin/env python
# coding: utf-8

# # Overview
# 
# This example demonstrates how to create and submit a job via [Job](https://docs.exabyte.io/api/Job/put_jobs_create) endpoints.

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
from settings import ENDPOINT_ARGS, ACCOUNT_ID
from utils import ensure_packages_are_installed, display_JSON
ensure_packages_are_installed()

from exabyte_api_client.endpoints.jobs import JobEndpoints
from exabyte_api_client.endpoints.materials import MaterialEndpoints
from exabyte_api_client.endpoints.workflows import WorkflowEndpoints


# ## Initialize the endpoints

# In[]:


job_endpoints = JobEndpoints(*ENDPOINT_ARGS)
material_endpoints = MaterialEndpoints(*ENDPOINT_ARGS)
workflow_endpoints = WorkflowEndpoints(*ENDPOINT_ARGS)


# Set job name.

# In[]:


JOB_NAME = "TEST JOB"


# ## Retrieve IDs
# 
# Default account's materail and workflow are used in this example to create the job. Adjust the queries to use different material and workflow.

# In[]:


default_material = material_endpoints.list({"isDefault": True, "owner._id": ACCOUNT_ID})[0]
default_workflow = workflow_endpoints.list({"isDefault": True, "owner._id": ACCOUNT_ID})[0]

material_id = default_material["_id"]
workflow_id = default_workflow["_id"]
owner_id = default_material["owner"]["_id"]


# ## Create job config
# 
# The job belongs to user's default account and it is created inside the defauult account's project. 

# In[]:


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

# In[]:


job = job_endpoints.create(config)
job_endpoints.submit(job['_id'])


# ## Print the job
# 
# Print the job in pretty JSON below. Check `status` field to make sure job is submiited.

# In[]:


job = job_endpoints.get(job['_id'])
display_JSON(job)

