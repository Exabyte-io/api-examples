#!/usr/bin/env python
# coding: utf-8

# # Overview
# 
# Inside this example we contact [Workflow](https://docs.exabyte.io/api/Workflows/get_workflows) endpoint to obtain a list of workflows that an account has access to.

# 1. Import required packages. Adjust [settings](../settings.ipynb) as necessary.

# In[1]:


import json
import argparse

from settings import *
from endpoints.workflows import WorkflowEndpoints


# 2. Set Parameters:
#     - **QUERY**: A query describing the documents to find. See [this](https://docs.meteor.com/api/collections.html#Mongo-Collection-find) for more information. 
#     - **OPTIONS**:
#         - **limit**: Maximum number of results to return. See [this](https://docs.meteor.com/api/collections.html#Mongo-Collection-find) for more information.

# In[9]:


QUERY = {
    "owner.slug": "demo"
}

OPTIONS = {
    "limit": 2
}


# 3. Initialize a helper class to interact with `WorkflowEndpoints`. This only has to be done once.

# In[10]:


endpoint = WorkflowEndpoints(HOST, PORT, ACCOUNT_ID, AUTH_TOKEN, VERSION, SECURE)


# 4. Contact the endpoint to list workflows according to the query above.

# In[11]:


workflows = endpoint.list(QUERY, OPTIONS)


# 5. Print the list of workflows saved under the corresponding variable in pretty JSON below.

# In[12]:


print json.dumps(workflows, indent=4)

