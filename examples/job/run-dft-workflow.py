#!/usr/bin/env python
# coding: utf-8

# # Overview
# 
# This example demonstrates how to use Exabyte RESTful API to study materials' properties by following the below steps:
# 
# - Importing materials from materials project
# 
# - Grouping imported materials inside a [materials set]()
# 
# - Creating jobs for the materials and grouping them inside a [jobs set]()
# 
# - Submiting jobs and monitoring the progress
# 
# - Extracting the [final structure]() and its properties
# 
# - Outputing the results as pandas dataFrame

# 1. Import required packages. Adjust [settings](../settings.ipynb) as necessary.

# In[24]:


import time
import pandas as pd
from IPython.display import IFrame

from settings import *
from endpoints.jobs import JobEndpoints
from endpoints.utils import flatten_material
from endpoints.materials import MaterialEndpoints
from endpoints.raw_properties import RawPropertiesEndpoints


# 2. Setup parameters
#     
#     - **MATERIALS_PROJECT_IDS**: A list of material IDs you would like to import.
# 
#     - **TAGS**: A list of tags you want to assign to imported materials.
# 
#     - **OWNER_ID**: ID of account that materials and jobs belongs to.
# 
#     - **PROJECT_ID**: The ID of project you would like to create the jobs in.
# 
#     - **JOBS_SET_NAME**: The name of the jobs set. Defaults to set.
# 
#     - **MATERIALS_SET_NAME**: The name of the materials set. Defaults to set.

# In[25]:


JOBS_SET_NAME = "set"
JOB_PREFIX = "phase-iii"
MATERIALS_SET_NAME = "set"
OWNER_ID = "knJyjbqxww7kt4GpA"
PROJECT_ID = "qqzbhRckbgTDNtD33"
TAGS = ["phase-iii", "difficulty-1"]
MATERIALS_PROJECT_IDS = ["mp-10694", "mp-29803"]


# This example is based on the below workflow. Make sure to copy it into your account and use its ID here. 

# In[5]:


IFrame('https://platform.exabyte.io/exabyte-io/workflows/YorcAEzG9gqtdvARj', width=800, height=650)


# In[26]:


WORKFLOW_ID = "5bdcdab3d7aad220b93140da"


# 3. Setup compute parameters
# 
#     - **NODES**: Number of nodes. Defaults to 1.
#     - **PPN**: Number of MPI processes per each node, Defaults to 1.
#     - **QUEUE**: The name of queue to submit the jobs into. Defaults to D.
#     - **TIME_LIMIT**: Job walltime. Defaults to "01:00:00" (one hour).
#     - **CLUSTER**: The full qualified domain name (FQDN) of the cluster to submit the jobs into:
#         - master-production-20160630-cluster-001.exabyte.io (AWS, c4, hyperthreaded)
#         - master-production-20160630-cluster-007.exabyte.io (Azure, H and NCv2 series, premium)

# In[64]:


PPN = "1"
QUEUE = "D"
NODES = "1"
TIME_LIMIT = "01:00:00"
CLUSTER = "master-production-20160630-cluster-001.exabyte.io"


# 4. Initialize job, material, and raw property endpoints. This needs to be done only once.

# In[28]:


job_endpoints = JobEndpoints(HOST, PORT, ACCOUNT_ID, AUTH_TOKEN, VERSION, SECURE)
material_endpoints = MaterialEndpoints(HOST, PORT, ACCOUNT_ID, AUTH_TOKEN, VERSION, SECURE)
raw_property_endpoints = RawPropertiesEndpoints(HOST, PORT, ACCOUNT_ID, AUTH_TOKEN, VERSION, SECURE)


# 5. Import material from materials project with the above tags.

# In[ ]:


materials = material_endpoints.import_from_materialsproject(MATERIALS_PROJECT_API_KEY, MATERIALS_PROJECT_IDS, OWNER_ID, TAGS)


# 6. Create a materials set and move the materials into it.

# In[31]:


materials_set = material_endpoints.create_set({"name": MATERIALS_SET_NAME, "owner": {"_id": OWNER_ID}})
for material in materials: material_endpoints.move_to_set(material["_id"], "", materials_set["_id"])


# 7. Create jobs for the materials above.

# In[32]:


compute = job_endpoints.get_compute(CLUSTER, PPN, NODES, QUEUE, TIME_LIMIT)
jobs = job_endpoints.create_by_ids(materials, WORKFLOW_ID, PROJECT_ID, OWNER_ID, JOB_PREFIX, compute)


# 8. Create a jobs set and move the jobs into it.

# In[33]:


jobs_set = job_endpoints.create_set({"name": JOBS_SET_NAME, "projectId": PROJECT_ID, "owner": {"_id": OWNER_ID}})
for job in jobs: job_endpoints.move_to_set(job["_id"], "", jobs_set["_id"])


# 9. Submit the jobs to the cluster.

# In[34]:


for job in jobs: job_endpoints.submit(job["_id"])


# 10. Wait for jobs to finish.

# In[35]:


job_ids = [job["_id"] for job in jobs]
while True:
    jobs = job_endpoints.list({"_id": {"$in": job_ids}}, {"status": 1})
    if all([job["status"] not in ["pre-submission", "submitted", "active"] for job in jobs]): break
    time.sleep(30)


# The following function returns a material property extacted in the given unit of the job's subworkflow. 

# In[36]:


def get_property_by_subworkow_and_unit_indecies(property_name, job, subworkflow_index, unit_index):
    unit_flowchart_id = job["workflow"]["subworkflows"][subworkflow_index]["units"][unit_index]["flowchartId"]
    return raw_property_endpoints.get_property(job["_id"], unit_flowchart_id, property_name)


# 11. For each material, extract final structure, pressure, and band gaps and flatten them to form the final dataFrame.

# In[39]:


rows = []
for material in materials:
    job = next((job for job in jobs if job["_material"]["_id"] == material["_id"]))
    final_structure = get_property_by_subworkow_and_unit_indecies("final_structure", job, 0, 0)["data"]
    pressure = get_property_by_subworkow_and_unit_indecies("pressure", job, 0, 0)["data"]["value"]
    unit_flowchart_id = job["workflow"]["subworkflows"][1]["units"][1]["flowchartId"]
    band_gap_direct = raw_property_endpoints.get_direct_band_gap(job["_id"], unit_flowchart_id)
    band_gap_indirect = raw_property_endpoints.get_indirect_band_gap(job["_id"], unit_flowchart_id)
    data = flatten_material(material)
    data.extend(flatten_material(final_structure))
    data.extend([pressure, band_gap_direct, band_gap_indirect])
    rows.append(data)


# 12. Create and print the final data as dataFrame.

# In[65]:


headers = []
keys = ["ID", "NAME", "TAGS", "N-SITES", "LAT-A", "LAT-B", "LAT-C", "LAT-ALPHA", "LAT-BETA", "LAT-GAMMA"]
headers.extend(["-".join(("INITIAL", key)) for key in keys])
headers.extend(["-".join(("FINAL", key)) for key in keys])
headers.extend(["PRESSURE", "DIRECT-GAP", "INDIRECT-GAP"])
pd.set_option('display.max_colwidth', -1)
pd.set_option('display.max_columns', 30)
pd.DataFrame(data=rows, columns=headers)

