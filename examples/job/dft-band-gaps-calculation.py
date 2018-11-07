#!/usr/bin/env python
# coding: utf-8

# # Overview
# 
# This example demonstrates how to use Exabyte RESTful API to create simulation [Jobs](https://docs.exabyte.io/jobs/overview/) programmatically for multiple [Materials](https://docs.exabyte.io/materials/overview/) at a time and extract the resulting [Properties](https://docs.exabyte.io/properties/overview/) forming a [Pandas](https://pandas.pydata.org/) dataframe.
# 
# This approach can work with any [Workflows](https://docs.exabyte.io/workflows/overview/) and multiple materials at once. For the demonstration purpose we use the Density Functional Theory and extract Electronic Band Gap as the property of interest.
# 
# 
# # Steps
# 
# We follow the below steps:
# 
# - Import materials from [materials project](https://materialsproject.org/)
# 
# - Group imported materials inside a [materials set](https://docs.exabyte.io/entities-general/sets/)
# 
# - Create jobs for the materials and grouping them inside a [jobs set](https://docs.exabyte.io/entities-general/sets/)
# 
# - Submit jobs and monitoring the progress
# 
# - Extract the final structure (relaxed structure) and its properties
# 
# - Output the results as Pandas dataFrame
# 
# # Pre-requisites
# 
# The explanation below assumes that the reader is familiar with the concepts used in Exabyte platform and RESTful API. We outline these below and direct the reader to the original sources of information:
# 
# - [Generating RESTful API authentication parameters](../system/get_authentication_params.ipynb)
# - [Importing materials from materials project](../material/import_materials_from_materialsproject.ipynb)
# - [Creating and submitting jobs](../job/create_and_submit_job.ipynb)

# # Execution
# 
# > <span style="color: orange">**NOTE**</span>: In order to run this example, an active Exabyte.io account with VASP access is required. RESTful API credentials shall be updated in [settings](../settings.ipynb). The generation of the credentials is also explained therein.
# 
# ## Import packages

# In[19]:


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
# Set the slug of [account](https://docs.exabyte.io/accounts/overview/) under which all the steps will be executed below.
# 
# > <span style="color: orange">**NOTE**</span>: The above step is required!

# In[20]:


ACCOUNT_SLUG = "exabyte"


# Set parameters for the materials to be imported:
#     
# - **MATERIALS_PROJECT_IDS**: a list of material IDs to be imported from materials project
# - **TAGS**: a list of [tags](https://docs.exabyte.io/entities-general/data/#tags) to assign to imported materials
# - **MATERIALS_SET_NAME**: the name of the materials set
# 

# In[21]:


MATERIALS_PROJECT_IDS = ["mp-10694", "mp-29803"]
MATERIALS_SET_NAME = "materials-set"
TAGS = ["tag1", "tag2"]


# Set parameters for the jobs to be ran for the imported materials:
# 
# - **JOB_NAME_PREFIX**: prefix to be used for the job name with "{JOB_NAME_PREFIX} {FORMULA}" convention (e.g.  "Job Name Prefix - SiGe")
# - **JOBS_SET_NAME**: the name of the jobs set
# - **PROJECT_SLUG**: slug of the [project](https://docs.exabyte.io/jobs/projects/) that the jobs will be created in. Below the default project ("Default") is used
# 

# In[22]:


PROJECT_SLUG = ACCOUNT_SLUG + "-default"
JOB_NAME_PREFIX = "Job Name Prefix"
JOBS_SET_NAME = "jobs-set"


# This example is based on [this](https://platform.exabyte.io/analytics/workflows/BEWfDREDFFL9g8Qpk) bank workflow which is later copied to the account workflows.

# In[23]:


BANK_WORKFLOW_ID = "BEWfDREDFFL9g8Qpk"


# Setup compute parameters. See [this](https://docs.exabyte.io/infrastructure/compute-settings/ui) for more information about compute parameters.
# 
# - **NODES**: Number of nodes. Defaults to 1.
# - **PPN**: Number of MPI processes per each node, Defaults to 1.
# - **QUEUE**: The name of queue to submit the jobs into. Defaults to D.
# - **TIME_LIMIT**: Job walltime. Defaults to "01:00:00" (one hour).
# - **CLUSTER**: The full qualified domain name (FQDN) of the cluster to submit the jobs into.

# In[33]:


PPN = "1"
QUEUE = "D"
NODES = "1"
TIME_LIMIT = "01:00:00"
CLUSTER = "master-production-20160630-cluster-001.exabyte.io"


# ## Initialize the endpoints

# In[25]:


job_endpoints = JobEndpoints(HOST, PORT, ACCOUNT_ID, AUTH_TOKEN, VERSION, SECURE)
project_endpoints = ProjectEndpoints(HOST, PORT, ACCOUNT_ID, AUTH_TOKEN, VERSION, SECURE)
material_endpoints = MaterialEndpoints(HOST, PORT, ACCOUNT_ID, AUTH_TOKEN, VERSION, SECURE)
raw_property_endpoints = RawPropertiesEndpoints(HOST, PORT, ACCOUNT_ID, AUTH_TOKEN, VERSION, SECURE)
bank_workflow_endpoints = BankWorkflowEndpoints(HOST, PORT, ACCOUNT_ID, AUTH_TOKEN, VERSION, SECURE)


# ## Create workflow
# 
# Retrieve account and project IDs as they are needed by the endpoints. 
# 
# Account's default material is used to extract the owner ID. You can extract the owner ID from any other account's [entities](https://docs.exabyte.io/entities-general/overview/).

# In[26]:


owner_id = material_endpoints.list({"isDefault": True, "owner.slug": ACCOUNT_SLUG})[0]["owner"]["_id"]
project_id = project_endpoints.list({"slug": PROJECT_SLUG, "owner.slug": ACCOUNT_SLUG})[0]["_id"]


# Copy bank workflow to the account's workflows.

# In[27]:


workflow_id = bank_workflow_endpoints.copy(BANK_WORKFLOW_ID, owner_id)["_id"]


# ## Import materials
# 
# Import materials from materials project with the above tags.

# In[28]:


materials = material_endpoints.import_from_materialsproject(MATERIALS_PROJECT_API_KEY, MATERIALS_PROJECT_IDS, owner_id, TAGS)


# Create a materials set and move the materials into it.

# In[29]:


materials_set = material_endpoints.create_set({"name": MATERIALS_SET_NAME, "owner": {"_id": owner_id}})
for material in materials: material_endpoints.move_to_set(material["_id"], "", materials_set["_id"])


# ## Create jobs
# 
# Create jobs for the materials above.

# In[34]:


compute = job_endpoints.get_compute(CLUSTER, PPN, NODES, QUEUE, TIME_LIMIT)
jobs = job_endpoints.create_by_ids(materials, workflow_id, project_id, owner_id, JOB_NAME_PREFIX, compute)


# Create a jobs set and move the jobs into it.

# In[35]:


jobs_set = job_endpoints.create_set({"name": JOBS_SET_NAME, "projectId": project_id, "owner": {"_id": owner_id}})
for job in jobs: job_endpoints.move_to_set(job["_id"], "", jobs_set["_id"])


# Submit the jobs for execution.

# In[36]:


for job in jobs: job_endpoints.submit(job["_id"])


# Monitor the jobs and print the status until they are all finished.

# In[37]:


job_ids = [job["_id"] for job in jobs]
wait_for_jobs_to_finish(job_endpoints, job_ids)


# ## Extract the results
# 
# For each material, extract final structure, pressure and band gaps. 
# 
# - Final structure and pressure are extracted from the first unit (vasp_relax with index 0) of the first job's subworkflow (volume-relaxation with index 0)
# 
# - Band gaps are extracted from the second unit (vasp-bands with index 1) of the second job's subworkflow (SCF-BS-BG-DOS with index 1).

# In[39]:


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
# 
# The below for-loop iterates over the results and flatten them to form the final Pandas dataFrame.

# In[41]:


table = []
for result in results:
    data = flatten_material(result["initial_structure"])
    data.extend(flatten_material(result["initial_structure"]))
    data.extend([result["pressure"], result["band_gap_direct"], result["band_gap_indirect"]])
    table.append(data)


# ## Ouput the results
# 
# Form the Pandas dataFrame headers according to the table generated above with the following abbreviations:
# 
# - **"INI"**: INITIAL
# - **"FIN"**: FINAL
# - **"N-SITES"**: Number of Sites
# - **"LAT"**: LATTICE

# In[42]:


headers = []
keys = ["ID", "NAME", "TAGS", "NS", "LAT-A", "LAT-B", "LAT-C", "LAT-ALPHA", "LAT-BETA", "LAT-GAMMA"]
headers.extend(["-".join(("INI", key)) for key in keys])
headers.extend(["-".join(("FIN", key)) for key in keys])
headers.extend(["PRESSURE", "DIRECT-GAP", "INDIRECT-GAP"])


# Create and print the final table as Pandas dataFrame.

# In[44]:


pd.set_option('display.max_colwidth', -1)
pd.set_option('display.max_columns', 30)
pd.DataFrame(data=table, columns=headers)

