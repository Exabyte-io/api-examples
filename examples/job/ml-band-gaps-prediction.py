#!/usr/bin/env python
# coding: utf-8

# # Overview
# 
# This example demonstrates how to use Exabyte RESTful API to build a machine learning (ML) model for a set of materials (named train materials) and use the model to predict the band gap of another set (named predict materials).
# 
# 
# 
# # Steps
# 
# We follow the below steps:
# 
# - Import materials from [materials project](https://materialsproject.org/)
# 
# - Calculate band gap for the train materials
# 
# - Build ML Train model for the train materials
# 
# - Create and submit a job to predict band gap for the predict materials
# 
# - Extract band gap for predict materials
# 
# - Output the results as Pandas dataFrame
# 
# # Pre-requisites
# 
# The explanation below assumes that the reader is familiar with the concepts used in Exabyte platform and RESTful API. We outline these below and direct the reader to the original sources of information:
# 
# - [Generating RESTful API authentication parameters](../system/get_authentication_params.ipynb)
# - [Importing materials from materials project](../material/import_materials_from_materialsproject.ipynb)
# - [Creating and submitting jobs](./create_and_submit_job.ipynb)
# - [Running DFT calculations](./dft-band-gaps-calculation.ipynb)

# # Execution
# 
# > <span style="color: orange">**NOTE**</span>: In order to run this example, an active Exabyte.io account is required. RESTful API credentials shall be updated in [settings](../settings.ipynb). The generation of the credentials is also explained therein.
# 
# ## Import packages

# In[52]:


import time
import pandas as pd
from IPython.display import IFrame

from utils import *
from settings import *
from endpoints.jobs import JobEndpoints
from endpoints.utils import flatten_material
from endpoints.projects import ProjectEndpoints
from endpoints.materials import MaterialEndpoints
from endpoints.workflows import WorkflowEndpoints
from endpoints.bank_workflows import BankWorkflowEndpoints
from endpoints.raw_properties import RawPropertiesEndpoints


# ## Setup parameters
# 
# Set the slug of [account](https://docs.exabyte.io/accounts/overview/) under which all the steps will be executed below.
# 
# > <span style="color: orange">**NOTE**</span>: The above step is required!

# In[53]:


ACCOUNT_SLUG = "exabyte"


# Set parameters for the materials to be imported:
#     
# - **TRAIN_MATERIALS_PROJECT_IDS**: a list of train material IDs to be imported from materials project
# - **PREDICT_MATERIALS_PROJECT_IDS**: a list of predict material IDs to be imported from materials project

# In[54]:


TRAIN_MATERIALS_PROJECT_IDS = ["mp-10694"]
PREDICT_MATERIALS_PROJECT_IDS = ["mp-29803"]


# Set parameters for the jobs to be ran for the imported materials:
# 
# - **JOB_NAME_PREFIX**: prefix to be used for the job name with "{JOB_NAME_PREFIX} {FORMULA}" convention (e.g.  "Job Name Prefix - SiGe")
# - **PROJECT_SLUG**: slug of the [project](https://docs.exabyte.io/jobs/projects/) that the jobs will be created in. Below the default project ("Default") is used
# 

# In[55]:


PROJECT_SLUG = ACCOUNT_SLUG + "-default"
JOB_NAME_PREFIX = "Job Name Prefix"


# Setup compute parameters. See [this](https://docs.exabyte.io/infrastructure/compute-settings/ui) for more information about compute parameters.
# 
# - **NODES**: Number of nodes. Defaults to 1.
# - **PPN**: Number of MPI processes per each node, Defaults to 1.
# - **QUEUE**: The name of queue to submit the jobs into. Defaults to D.
# - **TIME_LIMIT**: Job walltime. Defaults to "01:00:00" (one hour).
# - **CLUSTER**: The full qualified domain name (FQDN) or alias of the cluster to submit the jobs into.

# In[56]:


PPN = "1"
QUEUE = "D"
NODES = "1"
TIME_LIMIT = "01:00:00"
CLUSTER = "cluster-001"


# ## Initialize the endpoints

# In[57]:


job_endpoints = JobEndpoints(HOST, PORT, ACCOUNT_ID, AUTH_TOKEN, VERSION, SECURE)
project_endpoints = ProjectEndpoints(HOST, PORT, ACCOUNT_ID, AUTH_TOKEN, VERSION, SECURE)
material_endpoints = MaterialEndpoints(HOST, PORT, ACCOUNT_ID, AUTH_TOKEN, VERSION, SECURE)
workflow_endpoints = WorkflowEndpoints(HOST, PORT, ACCOUNT_ID, AUTH_TOKEN, VERSION, SECURE)
bank_workflow_endpoints = BankWorkflowEndpoints(HOST, PORT, ACCOUNT_ID, AUTH_TOKEN, VERSION, SECURE)
raw_property_endpoints = RawPropertiesEndpoints(HOST, PORT, ACCOUNT_ID, AUTH_TOKEN, VERSION, SECURE)


# ## Retrieve owner and project IDs
# 
# Retrieve owner and project IDs as they are needed by the endpoints. 
# 
# Account's default material is used to extract the owner ID. You can extract the owner ID from any other account's [entities](https://docs.exabyte.io/entities-general/overview/).

# In[58]:


owner_id = material_endpoints.list({"isDefault": True, "owner.slug": ACCOUNT_SLUG})[0]["owner"]["_id"]
project_id = project_endpoints.list({"slug": PROJECT_SLUG, "owner.slug": ACCOUNT_SLUG})[0]["_id"]


# ## Create workflows
# 
# Copy "ML: Train Model" and "Band Gap" bank workflows to the account's workflows.

# In[59]:


band_gap_workflow_id = copy_bank_workflow_by_system_name(bank_workflow_endpoints, "espresso-band-gap", owner_id)
ml_train_workflow_id = copy_bank_workflow_by_system_name(bank_workflow_endpoints, "exabyteml-ml-train-model", owner_id)


# ## Import materials
# 
# Import materials from materials project .

# In[61]:


train_materials = material_endpoints.import_from_materialsproject(MATERIALS_PROJECT_API_KEY, TRAIN_MATERIALS_PROJECT_IDS, owner_id)
predict_materials = material_endpoints.import_from_materialsproject(MATERIALS_PROJECT_API_KEY, PREDICT_MATERIALS_PROJECT_IDS, owner_id)


# ## Calculate band gap for train materials
# 
# Create jobs for the train materials.

# In[62]:


compute = job_endpoints.get_compute(CLUSTER, PPN, NODES, QUEUE, TIME_LIMIT)
jobs = job_endpoints.create_by_ids(train_materials, band_gap_workflow_id, project_id, owner_id, JOB_NAME_PREFIX, compute)


# Submit the jobs for execution.

# In[63]:


for job in jobs: job_endpoints.submit(job["_id"])


# Monitor the jobs and print the status until they are all finished.

# In[64]:


job_ids = [job["_id"] for job in jobs]
wait_for_jobs_to_finish(job_endpoints, job_ids)


# ## Build ML Train model
# 
# Create ML Train job for the train materials.

# In[65]:


name = "-".join((JOB_NAME_PREFIX, "train"))
material_ids = [m["_id"] for m in train_materials]
config = job_endpoints.get_config(material_ids, ml_train_workflow_id, project_id, owner_id, name, compute, True)
job = job_endpoints.create(config)


# Submit the train job for execution.

# In[66]:


job_endpoints.submit(job["_id"])


# Monitor the job and print the status until it is done.

# In[67]:


wait_for_jobs_to_finish(job_endpoints, [job["_id"]])


# ## Extract the predict workflow
# 
# Predict workflow is extracted from the last unit (train with index 4) of the first job's subworkflow (ML: Train Model with index 0).

# In[68]:


ml_predict_workflow_id = get_property_by_subworkow_and_unit_indecies(raw_property_endpoints, "workflow:ml_predict", job, 0, 4)["data"]["_id"]


# ## Create ML Predict job
# 
# Create ML Predict job for the predict materials.

# In[69]:


name = "-".join((JOB_NAME_PREFIX, "predict"))
material_ids = [m["_id"] for m in predict_materials]
config = job_endpoints.get_config(material_ids, ml_predict_workflow_id, project_id, owner_id, name, compute, True)
job = job_endpoints.create(config)


# Submit the predict job for execution.

# In[70]:


job_endpoints.submit(job["_id"])


# Monitor the job and print the status until its done.

# In[71]:


wait_for_jobs_to_finish(job_endpoints, [job["_id"]])


# ## Extract predicted properties
# 
# Predicted properties are extracted from the last unit (score with index 3) of the first job's subworkflow (ml_predict_subworkflow with index 0).

# In[72]:



predicted_properties = get_property_by_subworkow_and_unit_indecies(raw_property_endpoints, "predicted_properties", job, 0, 3)["data"]["values"]


# ## Flatten the results
# 
# The below for-loop iterates over the results and flatten them to form the final Pandas dataFrame.

# In[73]:


table = []
for exabyte_id, properties in predicted_properties.iteritems():
    material = next((m for m in predict_materials if m["exabyteId"] == exabyte_id))
    band_gaps = next((v for v in properties if v["name"] == "band_gaps"))
    direct_gap = next((v for v in band_gaps["values"] if v["type"] == "direct"))["value"]
    indirect_gap = next((v for v in band_gaps["values"] if v["type"] == "indirect"))["value"]
    table.append([material["_id"], material["name"], material["formula"], material["exabyteId"], direct_gap, indirect_gap])


# ## Ouput the results
# 
# Create and print the final table as Pandas dataFrame.

# In[74]:


headers = ["ID", "NAME", "FORMULA", "EXABYTE-ID", "DIRECT-GAP", "INDIRECT-GAP"]
pd.DataFrame(data=table, columns=headers)

