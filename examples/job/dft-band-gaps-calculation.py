#!/usr/bin/env python
# coding: utf-8

# # Overview
# 
# This example demonstrates how to use Exabyte RESTful API to create simulation Jobs programmatically for multiple Materials at a time and extract the resulting properties forming a Pandas dataframe.
# 
# This approach can work with any Workflow and multiple materials at once. For the demonstration purpose we use the Density Functional Theory and extract Electronic Band Gap as the property of interest.
# 
# 
# # Steps
# 
# We following the below steps:
# 
# - Import [materials](https://docs.exabyte.io/materials/overview/) from materials project
# 
# - Group imported materials inside a [materials set](https://docs.exabyte.io/entities-general/sets/)
# 
# - Create [jobs](https://docs.exabyte.io/jobs/overview/) for the materials and grouping them inside a [jobs set](https://docs.exabyte.io/entities-general/sets/)
# 
# - Submit jobs and monitoring the progress
# 
# - Extract the [final structure](https://docs.exabyte.io/properties/overview/) (relaxed structure) and its properties
# 
# - Output the results as [pandas](https://pandas.pydata.org/) dataFrame
# 
# # Pre-requisites
# 
# The explanation below assumes that the reader is familiar with the concepts used in Exabyte platform and RESTful API. We outline these below and direct the reader to the original sources of information:
# 
# - [Generating RESTFul API authentication parameters](../system/get_authentication_params.ipynb)
# - [Importing materials from materials project](../material/import_materials_from_materialsproject.ipynb)
# - [Creating and submitting jobs](../job/create_and_submit_job.ipynb)

# # Execution
# 
# > <span style="color: orange">**NOTE**</span>: In order to run this example, an active account with Exabyte.io is required. RESTful API credentials shall be updated in [settings](../settings.ipynb). The generation of the credentials is also explained therein.
# 
# ## Import packages

# In[1]:


import time
import pandas as pd
from IPython.display import IFrame

from utils import *
from settings import *
from endpoints.jobs import JobEndpoints
from endpoints.utils import flatten_material
from endpoints.projects import ProjectEndpoints
from endpoints.materials import MaterialEndpoints
from endpoints.bank_workflows import BankWorkflowEndpoints
from endpoints.raw_properties import RawPropertiesEndpoints


# ## Setup parameters
# 
# Set the account under which all the steps will be executed below:
# 
# - **ACCOUNT_SLUG**: Slug of [account](https://docs.exabyte.io/accounts/overview/) entities belong to.
# 
# > <span style="color: orange">**NOTE**</span>: The above step is required!

# In[18]:


ACCOUNT_SLUG = "exabyte"


# Set parameters for the materials to be imported:
#     
# - **MATERIALS_PROJECT_IDS**: A list of material IDs to be imported from [materials project](https://materialsproject.org/).
# - **TAGS**: A list of [tags](https://docs.exabyte.io/entities-general/actions/metadata/) to assign to imported materials.
# - **MATERIALS_SET_NAME**: The name of the materials set. Defaults to "materials-set".
# 

# In[2]:


MATERIALS_PROJECT_IDS = ["mp-10694", "mp-29803"]
MATERIALS_SET_NAME = "materials-set"
TAGS = ["phase-iii", "difficulty-1"]


# Set parameters for the jobs to be ran for the imported materials:
# 
# - **JOB_NAME_PREFIX**: prefix to be used for the jobs created for imported materials, the resulting names will be similar to "{JOB_NAME_PREFIX} {FORMULA}" - "Job Name Prefix - SiGe",
# - **JOBS_SET_NAME**: the name of the jobs set. Defaults to "jobs-set",
# - **PROJECT_SLUG**: slug of the [project](https://docs.exabyte.io/jobs/projects/) that the jobs will be created in. Below the default project ("Default") is used.
# 

# In[2]:


PROJECT_SLUG = ACCOUNT_SLUG + "-default"
JOB_NAME_PREFIX = "Job Name Prefix"
JOBS_SET_NAME = "jobs-set"


# This example is based on the below [bank workflow](https://platform.exabyte.io/analytics/workflows/BEWfDREDFFL9g8Qpk) which is later copied to the account workflows.

# In[3]:


BANK_WORKFLOW_ID = "BEWfDREDFFL9g8Qpk"


# Setup compute parameters:
# 
# - **NODES**: Number of nodes. Defaults to 1.
# - **PPN**: Number of MPI processes per each node, Defaults to 1.
# - **QUEUE**: The name of queue to submit the jobs into. Defaults to D.
# - **TIME_LIMIT**: Job walltime. Defaults to "01:00:00" (one hour).
# - **CLUSTER**: The full qualified domain name (FQDN) of the cluster to submit the jobs into:
#     - master-production-20160630-cluster-001.exabyte.io (AWS, c4, hyperthreaded)
#     - master-production-20160630-cluster-007.exabyte.io (Azure, H and NCv2 series, premium)

# In[4]:


PPN = "1"
QUEUE = "D"
NODES = "1"
TIME_LIMIT = "01:00:00"
CLUSTER = "master-production-20160630-cluster-001.exabyte.io"


# ## Initialize the endpoints

# In[5]:


job_endpoints = JobEndpoints(HOST, PORT, ACCOUNT_ID, AUTH_TOKEN, VERSION, SECURE)
project_endpoints = ProjectEndpoints(HOST, PORT, ACCOUNT_ID, AUTH_TOKEN, VERSION, SECURE)
material_endpoints = MaterialEndpoints(HOST, PORT, ACCOUNT_ID, AUTH_TOKEN, VERSION, SECURE)
raw_property_endpoints = RawPropertiesEndpoints(HOST, PORT, ACCOUNT_ID, AUTH_TOKEN, VERSION, SECURE)
bank_workflow_endpoints = BankWorkflowEndpoints(HOST, PORT, ACCOUNT_ID, AUTH_TOKEN, VERSION, SECURE)


# ## Create workflow
# 
# Obtain account and project IDs as the endpoints work with IDs than slugs.

# In[6]:


owner_id = material_endpoints.list({"isDefault": True, "owner.slug": ACCOUNT_SLUG})[0]["owner"]["_id"]
project_id = project_endpoints.list({"slug": PROJECT_SLUG, "owner.slug": ACCOUNT_SLUG})[0]["_id"]


# Copy bank workflow to account workflows.

# In[7]:


workflow_id = bank_workflow_endpoints.copy(BANK_WORKFLOW_ID, owner_id)["_id"]


# ## Import materials
# 
# Import material from materials project with the above tags.

# In[8]:


materials = material_endpoints.import_from_materialsproject(MATERIALS_PROJECT_API_KEY, MATERIALS_PROJECT_IDS, owner_id, TAGS)


# Create a materials set and move the materials into it.

# In[9]:


materials_set = material_endpoints.create_set({"name": MATERIALS_SET_NAME, "owner": {"_id": owner_id}})
for material in materials: material_endpoints.move_to_set(material["_id"], "", materials_set["_id"])


# ## Create jobs
# 
# Create jobs for the materials above.

# In[10]:


compute = job_endpoints.get_compute(CLUSTER, PPN, NODES, QUEUE, TIME_LIMIT)
jobs = job_endpoints.create_by_ids(materials, workflow_id, project_id, owner_id, JOB_PREFIX, compute)


# Create a jobs set and move the jobs into it.

# In[11]:


jobs_set = job_endpoints.create_set({"name": JOBS_SET_NAME, "projectId": project_id, "owner": {"_id": owner_id}})
for job in jobs: job_endpoints.move_to_set(job["_id"], "", jobs_set["_id"])


# Submit the jobs for execution

# In[12]:


for job in jobs: job_endpoints.submit(job["_id"])


# Wait for jobs to finish.

# In[13]:


job_ids = [job["_id"] for job in jobs]
wait_for_jobs_to_finish(job_endpoints, job_ids)


# The following function returns a material property extracted in the given unit of the job's subworkflow. 

# In[14]:


def get_property_by_subworkow_and_unit_indecies(property_name, job, subworkflow_index, unit_index):
    """
    Returns the property extracted in the given unit of the job's subworkflow. 
    """
    unit_flowchart_id = job["workflow"]["subworkflows"][subworkflow_index]["units"][unit_index]["flowchartId"]
    return raw_property_endpoints.get_property(job["_id"], unit_flowchart_id, property_name)


# ## Extract the results
# 
# For each material, extract final structure, pressure and add them to "results".

# In[15]:


results = []
for material in materials:
    job = next((job for job in jobs if job["_material"]["_id"] == material["_id"]))
    final_structure = get_property_by_subworkow_and_unit_indecies("final_structure", job, 0, 0)["data"]
    pressure = get_property_by_subworkow_and_unit_indecies("pressure", job, 0, 0)["data"]["value"]
    unit_flowchart_id = job["workflow"]["subworkflows"][1]["units"][1]["flowchartId"]
    band_gap_direct = raw_property_endpoints.get_direct_band_gap(job["_id"], unit_flowchart_id)
    band_gap_indirect = raw_property_endpoints.get_indirect_band_gap(job["_id"], unit_flowchart_id)
    results.append({
        "initial_structure": material,
        "final_structure": final_structure,
        "pressure": pressure,
        "band_gap_direct": band_gap_direct,
        "band_gap_indirect": band_gap_indirect,
    })


# ## Flatten the results

# In[ ]:


table = []
for result in results:
    data = flatten_material(result["initial_structure"])
    data.extend(flatten_material(result["initial_structure"]))
    data.extend(result["pressure"], result["band_gap_direct"], result["band_gap_indirect"])
    table.append(data)


# ## Ouput the results
# 
# Form the headers with the following abbreviations:
# 
# - **"INI"**: INITIAL
# - **"FIN"**: FINAL
# - **"NS"**: Number of sites
# - **"LAT"**: LATTICE

# In[ ]:


headers = []
keys = ["ID", "NAME", "TAGS", "NS", "LAT-A", "LAT-B", "LAT-C", "LAT-ALPHA", "LAT-BETA", "LAT-GAMMA"]
headers.extend(["-".join(("INI", key)) for key in keys])
headers.extend(["-".join(("FIN", key)) for key in keys])
headers.extend(["PRESSURE", "DIRECT-GAP", "INDIRECT-GAP"])


# Create and print the final data as dataFrame.

# In[16]:


headers = []
data = [[flatten_material(d["initial_structure"]), flatten_material(d["initial_structure"], d["pressure"], d["band_gap_direct"], d["band_gap_indirect"]] for d in data]
keys = ["ID", "NAME", "TAGS", "N-SITES", "LAT-A", "LAT-B", "LAT-C", "LAT-ALPHA", "LAT-BETA", "LAT-GAMMA"]
headers.extend(["-".join(("INITIAL", key)) for key in keys])
headers.extend(["-".join(("FINAL", key)) for key in keys])
headers.extend(["PRESSURE", "DIRECT-GAP", "INDIRECT-GAP"])
pd.set_option('display.max_colwidth', -1)
pd.set_option('display.max_columns', 30)
pd.DataFrame(data=table, columns=headers)

