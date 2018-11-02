#!/usr/bin/env python
# coding: utf-8

# # Overview
# 
# This example demonstrates how to create and submit a job via [Job](https://docs.exabyte.io/api/Job/put_jobs_create) endpoints.

# 1. Import required packages. Adjust [settings](../settings.ipynb) as necessary.

# In[2]:


import json
import argparse
import nbimporter
nbimporter.options['only_defs'] = False

from settings import *
from endpoints.jobs import JobEndpoints


# 2. Create job config in JSON format. Adjust material ID, workflow ID and job name accordingly.

# In[4]:


CONFIG = {
    "_material": {
        "_id": "tkmEX2KdrSzgo54km"
    },
    "workflow": {
        "_id": "6bd94da34794abf42a697fe1"
    },
    "name": "TEST JOB"
}


# 3. Initialize `JobEndpoints` class and call `create` function to create the job.

# In[6]:


endpoint = JobEndpoints(HOST, PORT, ACCOUNT_ID, AUTH_TOKEN, VERSION, SECURE)
job = endpoint.create(CONFIG)


# 4. Submit job by its ID.

# In[7]:


endpoint.submit(job['_id'])


# 4. Print the job in pretty JSON below. Check `status` field to make sure job is submiited. 

# In[9]:


print json.dumps(endpoint.get(job['_id']), indent=4)

