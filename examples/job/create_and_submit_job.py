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
from exabyte_api_client.endpoints.jobs import JobEndpoints
from exabyte_api_client.endpoints.materials import MaterialEndpoints
from exabyte_api_client.endpoints.workflows import WorkflowEndpoints

# Import settings file
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path: sys.path.append(module_path)
from settings import ENDPOINT_ARGS, ACCOUNT_ID


# ## Initialize the endpoints

# In[15]:


job_endpoints = JobEndpoints(*ENDPOINT_ARGS)
material_endpoints = MaterialEndpoints(*ENDPOINT_ARGS)
workflow_endpoints = WorkflowEndpoints(*ENDPOINT_ARGS)


# Set job name.

# In[17]:


JOB_NAME = "TEST JOB"


# ## Retrieve IDs
# 
# Default account's materail and workflow are used in this example to create the job. Adjust the queries to use different material and workflow.

# In[18]:


default_material = material_endpoints.list({"isDefault": True, "owner._id": ACCOUNT_ID})[0]
default_workflow = workflow_endpoints.list({"isDefault": True, "owner._id": ACCOUNT_ID})[0]

material_id = default_material["_id"]
workflow_id = default_workflow["_id"]
owner_id = default_material["owner"]["_id"]


# ## Create job config
# 
# The job belongs to user's default account and it is created inside the defauult account's project. 

# In[19]:


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

# In[20]:


job = job_endpoints.create(config)
job_endpoints.submit(job['_id'])


# ## Print the job
# 
# Print the job in pretty JSON below. Check `status` field to make sure job is submiited.

# In[22]:


job = job_endpoints.get(job['_id'])
print(json.dumps(job, indent=4))

