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

# In[2]:


from endpoints.jobs import JobEndpoints
from endpoints.materials import MaterialEndpoints
from endpoints.workflows import WorkflowEndpoints
from settings import HOST, PORT, ACCOUNT_ID, AUTH_TOKEN, VERSION, SECURE


# ## Initialize the endpoints

# In[6]:


args = [HOST, PORT, ACCOUNT_ID, AUTH_TOKEN, VERSION, SECURE]
job_endpoints = JobEndpoints(*args)
material_endpoints = MaterialEndpoints(*args)
workflow_endpoints = WorkflowEndpoints(*args)


# ## Setup parameters
# 
# Set the slug of [account](https://docs.exabyte.io/accounts/overview/) under which all the steps will be executed below.
# 
# > <span style="color: orange">**NOTE**</span>: This step is mandatory!

# In[ ]:


ACCOUNT_SLUG = "exabyte"


# Set job name.

# In[ ]:


JOB_NAME = "TEST JOB"


# ## Retrieve material and workflow IDs
# 
# Default account's materail and workflow are used in this example to create the job. Adjust the queries to use different material and workflow.

# In[ ]:


material_id = material_endpoints.list({"isDefault": True, "owner.slug": ACCOUNT_SLUG})[0]["_id"]
workflow_id = material_endpoints.list({"isDefault": True, "owner.slug": ACCOUNT_SLUG})[0]["_id"]


# ## Create job config
# 
# The job belongs to user's default account and it is created inside the defauult account's project. 

# In[4]:


config = {
    "_material": {
        "_id": material_id
    },
    "workflow": {
        "_id": workflow_id
    },
    "name": JOB_NAME
}


# ## Create and submit job

# In[7]:


job = endpoint.create(config)
endpoint.submit(job['_id'])


# ## Print the job
# 
# Print the job in pretty JSON below. Check `status` field to make sure job is submiited.

# In[9]:


job = endpoint.get(job['_id'])
print json.dumps(job, indent=4)

