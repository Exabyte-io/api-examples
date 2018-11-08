#!/usr/bin/env python
# coding: utf-8

# # Overview
# 
# This example demonstrates how to create and submit a job via [Job](https://docs.exabyte.io/api/Job/put_jobs_create) endpoints.

# # Execution
# 
# > <span style="color: orange">**NOTE**</span>: In order to run this example, an active Exabyte.io account is required. RESTful API credentials shall be updated in [settings](../settings.ipynb). The generation of the credentials is also explained therein.
# 
# ## Import packages

# In[14]:


import json

from endpoints.jobs import JobEndpoints
from endpoints.materials import MaterialEndpoints
from endpoints.workflows import WorkflowEndpoints
from settings import HOST, PORT, ACCOUNT_ID, AUTH_TOKEN, VERSION, SECURE


# ## Initialize the endpoints

# In[15]:


args = [HOST, PORT, ACCOUNT_ID, AUTH_TOKEN, VERSION, SECURE]
job_endpoints = JobEndpoints(*args)
material_endpoints = MaterialEndpoints(*args)
workflow_endpoints = WorkflowEndpoints(*args)


# ## Setup parameters
# 
# Set the slug of [account](https://docs.exabyte.io/accounts/overview/) under which all the steps will be executed below.
# 
# > <span style="color: orange">**NOTE**</span>: This step is mandatory!

# In[16]:


ACCOUNT_SLUG = "exabyte"


# Set job name.

# In[17]:


JOB_NAME = "TEST JOB"


# ## Retrieve IDs
# 
# Default account's materail and workflow are used in this example to create the job. Adjust the queries to use different material and workflow.

# In[18]:


default_material = material_endpoints.list({"isDefault": True, "owner.slug": ACCOUNT_SLUG})[0]
default_workflow = workflow_endpoints.list({"isDefault": True, "owner.slug": ACCOUNT_SLUG})[0]

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
print json.dumps(job, indent=4)

