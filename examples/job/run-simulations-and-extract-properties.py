#!/usr/bin/env python
# coding: utf-8

# # Run Simulations and Extract Properties
# 
# This example demonstrates how to use Exabyte RESTful API to create simulation [Jobs](https://docs.exabyte.io/jobs/overview/) programmatically for multiple [Materials](https://docs.exabyte.io/materials/overview/) at once and extract the resulting [Properties](https://docs.exabyte.io/properties/overview/) forming a [Pandas](https://pandas.pydata.org/) dataframe.
# 
# This approach can work with any [Workflows](https://docs.exabyte.io/workflows/overview/). For the demonstration purpose we use the Density Functional Theory and extract Electronic Band Gap as the property of interest.
# 
# > <span style="color: orange">**IMPORTANT NOTE**</span>: In order to run this example in full, an active Exabyte.io account with access to VASP (Vienna ab-initio simulations package) is required. Alternatively, Readers may substitute the workflow ID below with another one (an equivalent one for Quantum ESPRESSO, for example) and adjust extraction of the results ("Extract results" section). RESTful API credentials shall be updated in [settings](../settings.py).
# 
# 
# ## Steps
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
# - Extract the [final structure](https://docs.exabyte.io/properties/structural/final-structure) (relaxed structure) and its properties
# 
# - Output the results as Pandas dataFrame
# 
# ## Pre-requisites
# 
# The explanation below assumes that the reader is familiar with the concepts used in Exabyte platform and RESTful API. We outline these below and direct the reader to the original sources of information:
# 
# - [Generating RESTful API authentication parameters](../system/get_authentication_params.ipynb)
# - [Importing materials from materials project](../material/import_materials_from_materialsproject.ipynb)
# - [Creating and submitting jobs](../job/create_and_submit_job.ipynb)

# ## Execution
# 
# 
# ### Import packages

# In[]:


import time
from IPython.display import IFrame
import os
import sys

# Import settings file and utils file
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path: sys.path.append(module_path)
from settings import ENDPOINT_ARGS, ACCOUNT_ID, MATERIALS_PROJECT_API_KEY
from utils import wait_for_jobs_to_finish, get_property_by_subworkow_and_unit_indicies, dataframe_to_html, ensure_packages_are_installed
ensure_packages_are_installed()

import pandas as pd

# Relevant functions from the API client
from exabyte_api_client.endpoints.jobs import JobEndpoints
from exabyte_api_client.utils.materials import flatten_material
from exabyte_api_client.endpoints.projects import ProjectEndpoints
from exabyte_api_client.endpoints.materials import MaterialEndpoints
from exabyte_api_client.endpoints.bank_workflows import BankWorkflowEndpoints
from exabyte_api_client.endpoints.raw_properties import RawPropertiesEndpoints


# #### Materials
#     
# - **MATERIALS_PROJECT_IDS**: a list of material IDs to be imported from materials project
# - **TAGS**: a list of [tags](https://docs.exabyte.io/entities-general/data/#tags) to assign to imported materials
# - **MATERIALS_SET_NAME**: the name of the materials set
# 

# In[]:


MATERIALS_PROJECT_IDS = ["mp-149", "mp-32"] # Si and Ge
MATERIALS_SET_NAME = "materials-set"
TAGS = ["tag1", "tag2"]


# #### Jobs
# 
# Parameters for the jobs to be ran for the imported materials:
# 
# - **JOB_NAME_PREFIX**: prefix to be used for the job name with "{JOB_NAME_PREFIX} {FORMULA}" convention (e.g.  "Job Name Prefix - SiGe")
# - **JOBS_SET_NAME**: the name of the jobs set

# In[]:


JOB_NAME_PREFIX = "Job Name Prefix"
JOBS_SET_NAME = "jobs-set"


# #### Workflow
# 
# This example is based on [this](https://platform.exabyte.io/analytics/workflows/56xDyXsPMNJ7cF9nv) bank workflow which is later copied to the account workflows collection.  The workflow is named "D3-GGA-BS-BG-DOS-ALL" and utilizes the logic explained in https://arxiv.org/pdf/1808.05325.pdf, for example (see section "Methodology", Table I). "D3" indicates the difficulty level 3 per the table convention. BS, BG, DOS indicate the properties extracted - Band Structure, Band Gap, Density of States. The workflow is utilizing VASP simulation engine at version 5.4.4.

# In[]:


BANK_WORKFLOW_ID = "56xDyXsPMNJ7cF9nv"


# In[]:


# Visualize the bank workflow below
# NOTE: might not be rendered on Github
IFrame("https://platform.exabyte.io/analytics/workflows/{}".format(BANK_WORKFLOW_ID), width=900, height=650)


# #### Compute
# 
# Setup compute parameters. See [this](https://docs.exabyte.io/infrastructure/compute-settings/ui) for more information about compute parameters.
# 
# - **NODES**: Number of nodes. Defaults to 1.
# - **PPN**: Number of MPI processes per each node, Defaults to 1.
# - **QUEUE**: The name of queue to submit the jobs into. Defaults to D.
# - **TIME_LIMIT**: Job walltime. Defaults to "01:00:00" (one hour).
# - **CLUSTER**: The full qualified domain name (FQDN) or alias of the cluster to submit the jobs into.

# In[]:


PPN = "1"
QUEUE = "D"
NODES = "1"
TIME_LIMIT = "01:00:00"
CLUSTER = "cluster-001"


# ### Initialize endpoints

# In[]:


job_endpoints = JobEndpoints(*ENDPOINT_ARGS)
project_endpoints = ProjectEndpoints(*ENDPOINT_ARGS)
material_endpoints = MaterialEndpoints(*ENDPOINT_ARGS)
raw_property_endpoints = RawPropertiesEndpoints(*ENDPOINT_ARGS)
bank_workflow_endpoints = BankWorkflowEndpoints(*ENDPOINT_ARGS)


# Next, we retrieve the owner and project IDs as they are needed by the endpoints. Account's default material is used to extract the owner ID. One can extract the owner ID from any other account's [entities](https://docs.exabyte.io/entities-general/overview/).

# In[]:


owner_id = material_endpoints.list({"isDefault": True, "owner._id": ACCOUNT_ID})[0]["owner"]["_id"]
project_id = project_endpoints.list({"isDefault": True, "owner._id": ACCOUNT_ID})[0]["_id"]


# ### Create workflow
# 
# Copy bank workflow (template) to the account's workflows collection.

# In[]:


workflow_id = bank_workflow_endpoints.copy(BANK_WORKFLOW_ID, owner_id)["_id"]


# ### Import materials
# 
# Import materials from materials project with the above tags.

# In[]:


materials = material_endpoints.import_from_materialsproject(MATERIALS_PROJECT_API_KEY, MATERIALS_PROJECT_IDS, owner_id, TAGS)


# Create a materials set and move the materials into it.

# In[]:


materials_set = material_endpoints.create_set({"name": MATERIALS_SET_NAME, "owner": {"_id": owner_id}})
for material in materials: material_endpoints.move_to_set(material["_id"], "", materials_set["_id"])


# ### Create jobs
# 
# Create jobs for the materials above.

# In[]:


compute = job_endpoints.get_compute(CLUSTER, PPN, NODES, QUEUE, TIME_LIMIT)
jobs = job_endpoints.create_by_ids(materials, workflow_id, project_id, owner_id, JOB_NAME_PREFIX, compute)


# Create a jobs set and move the jobs into it.

# In[]:


jobs_set = job_endpoints.create_set({"name": JOBS_SET_NAME, "projectId": project_id, "owner": {"_id": owner_id}})
for job in jobs: job_endpoints.move_to_set(job["_id"], "", jobs_set["_id"])


# Submit the jobs for execution.

# In[]:


for job in jobs: job_endpoints.submit(job["_id"])


# Monitor the jobs and print the status until they are all finished.

# In[]:


job_ids = [job["_id"] for job in jobs]
wait_for_jobs_to_finish(job_endpoints, job_ids)


# ### Extract results
# 
# For each material, simulaion job, final structure, pressure and band gaps are extracted. 
# 
# - Final structure and pressure are extracted from the first unit (vasp_relax with index 0) of the first job's subworkflow (volume-relaxation with index 0)
# 
# - Band gaps are extracted from the second unit (vasp-bands with index 1) of the second job's subworkflow (SCF-BS-BG-DOS with index 1).

# In[]:


results = []
for material in materials:
    job = next((job for job in jobs if job["_material"]["_id"] == material["_id"]))
    final_structure = get_property_by_subworkow_and_unit_indicies(raw_property_endpoints, "final_structure", job, 0, 0)["data"]
    pressure = get_property_by_subworkow_and_unit_indicies(raw_property_endpoints, "pressure", job, 0, 0)["data"]["value"]
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


# ### Flatten results
# 
# The below for-loop iterates over the results and flatten them to form the final Pandas dataFrame.

# In[]:


table = []
for result in results:
    data = flatten_material(result["initial_structure"])
    data.extend(flatten_material(result["initial_structure"]))
    data.extend([result["pressure"], result["band_gap_direct"], result["band_gap_indirect"]])
    table.append(data)


# ### Output results
# 
# Form the Pandas dataFrame headers according to the table generated above with the following abbreviations:
# 
# - **"INI"**: INITIAL
# - **"FIN"**: FINAL
# - **"N-SITES"**: Number of Sites
# - **"LAT"**: LATTICE

# In[]:


headers = []
keys = ["ID", "NAME", "TAGS", "NS", "LAT-A", "LAT-B", "LAT-C", "LAT-ALPHA", "LAT-BETA", "LAT-GAMMA"]
headers.extend(["-".join(("INI", key)) for key in keys])
headers.extend(["-".join(("FIN", key)) for key in keys])
headers.extend(["PRESSURE", "DIRECT-GAP", "INDIRECT-GAP"])


# Create and print the final table as Pandas dataFrame.

# In[]:


df = pd.DataFrame(data=table, columns=headers)
html = dataframe_to_html(df)
html

