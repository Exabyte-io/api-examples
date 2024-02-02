#!/usr/bin/env python
# coding: utf-8

# <a href="https://colab.research.google.com/github/Exabyte-io/api-examples/blob/dev/examples/job/ml-train-model-predict-properties.ipynb" target="_parent">
# <img alt="Open in Google Colab" src="https://user-images.githubusercontent.com/20477508/128780728-491fea90-9b23-495f-a091-11681150db37.jpeg" width="150" border="0">
# </a>

# # Overview
#
# This example demonstrates how to use Mat3ra RESTful API to build a machine learning (ML) model for a set of materials called "train materials" and use the model to predict properties of another set called "target materials". The general approach can work for multiple properties, we use the Electronic Band Gap in this example.
#
#
#
# ## Steps
#
# We follow the below steps:
#
# - Import materials from [materials project](https://materialsproject.org/)
# - Calculate band gap for the "train materials"
# - Build ML Train model based on the "train materials"
# - Create and submit a job to predict band gap for the "target materials"
# - Extract band gap for "target materials"
# - Output the results as Pandas dataFrame
#
# ## Pre-requisites
#
# The explanation below assumes that the reader is familiar with the concepts used in Mat3ra platform and RESTful API. We outline these below and direct the reader to the original sources of information:
#
# - [Generating RESTful API authentication parameters](../system/get_authentication_params.ipynb)
# - [Importing materials from materials project](../material/import_materials_from_materialsproject.ipynb)
# - [Creating and submitting jobs](./create_and_submit_job.ipynb)
# - [Running DFT calculations](./run-simulations-and-extract-properties.ipynb)

# # Complete Authorization Form and Initialize Settings
#
# This will also determine environment and set all environment variables. We determine if we are using Jupyter Notebooks or Google Colab to run this tutorial.
#
# If you are running this notebook from Google Colab, Colab takes ~1 min to execute the following cell.
#
# ACCOUNT_ID and AUTH_TOKEN - Authentication parameters needed for when making requests to [Mat3ra.com's API Endpoints](https://docs.mat3ra.com/rest-api/endpoints/).
#
# MATERIALS_PROJECT_API_KEY - Authentication parameter needed for when making requests to [Material Project's API](https://materialsproject.org/open)
#
# ORGANIZATION_ID - Authentication parameter needed for when working with collaborative accounts https://docs.mat3ra.com/collaboration/organizations/overview/
#
# > <span style="color: orange">**NOTE**</span>: If you are running this notebook from Jupyter, the variables ACCOUNT_ID, AUTH_TOKEN, MATERIALS_PROJECT_API_KEY, and ORGANIZATION_ID should be set in the file [settings.json](../../utils/settings.json) if you need to use these variables. To obtain API token parameters, please see the following link to the documentation explaining how to get them: https://docs.mat3ra.com/accounts/ui/preferences/api/

# In[ ]:


# @title Authorization Form
ACCOUNT_ID = "ACCOUNT_ID"  # @param {type:"string"}
AUTH_TOKEN = "AUTH_TOKEN"  # @param {type:"string"}
MATERIALS_PROJECT_API_KEY = "MATERIALS_PROJECT_API_KEY"  # @param {type:"string"}
ORGANIZATION_ID = "ORGANIZATION_ID"  # @param {type:"string"}

import os

if "COLAB_JUPYTER_IP" in os.environ:
    os.environ.update(
        dict(
            ACCOUNT_ID=ACCOUNT_ID,
            AUTH_TOKEN=AUTH_TOKEN,
            MATERIALS_PROJECT_API_KEY=MATERIALS_PROJECT_API_KEY,
            ORGANIZATION_ID=ORGANIZATION_ID,
        )
    )

    get_ipython().system('GIT_BRANCH="dev"; export GIT_BRANCH; curl -s "https://raw.githubusercontent.com/Exabyte-io/api-examples/${GIT_BRANCH}/scripts/env.sh" | bash')


# # Imports

# In[ ]:


import time
from IPython.display import IFrame

# Import settings file and utils file
from utils.settings import ENDPOINT_ARGS, ACCOUNT_ID, MATERIALS_PROJECT_API_KEY
from utils.generic import (
    dataframe_to_html,
    copy_bank_workflow_by_system_name,
    wait_for_jobs_to_finish,
    get_property_by_subworkflow_and_unit_indicies,
    display_JSON,
)

import pandas as pd

# Import relevant portions of the API client
from exabyte_api_client.endpoints.jobs import JobEndpoints
from exabyte_api_client.utils.materials import flatten_material
from exabyte_api_client.endpoints.projects import ProjectEndpoints
from exabyte_api_client.endpoints.materials import MaterialEndpoints
from exabyte_api_client.endpoints.workflows import WorkflowEndpoints
from exabyte_api_client.endpoints.bank_workflows import BankWorkflowEndpoints
from exabyte_api_client.endpoints.properties import PropertiesEndpoints


# #### Materials
#
# Set parameters for the materials to be imported:
#
# - **TRAIN_MATERIALS_PROJECT_IDS**: a list of material IDs to train ML model based on
# - **TARGET_MATERIALS_PROJECT_IDS**: a list of material IDs to predict the property for

# In[ ]:


TRAIN_MATERIALS_PROJECT_IDS = ["mp-149", "mp-978534"]  # Si, SiGe
TARGET_MATERIALS_PROJECT_IDS = ["mp-32"]  # Ge


# #### Jobs
#
# Set parameters for the jobs to be ran for the imported materials:
#
# - **JOB_NAME_PREFIX**: prefix to be used for the job name with "{JOB_NAME_PREFIX} {FORMULA}" convention (e.g.  "Job Name Prefix - SiGe")

# In[ ]:


JOB_NAME_PREFIX = "Job Name Prefix"


# #### Compute
#
# Setup compute parameters. See [this](https://docs.mat3ra.com/infrastructure/compute-settings/ui) for more information about compute parameters.
#
# - **NODES**: Number of nodes. Defaults to 1.
# - **PPN**: Number of MPI processes per each node, Defaults to 1.
# - **QUEUE**: The name of queue to submit the jobs into. Defaults to D.
# - **TIME_LIMIT**: Job walltime. Defaults to "01:00:00" (one hour).
# - **CLUSTER**: The full qualified domain name (FQDN) or alias of the cluster to submit the jobs into.

# In[ ]:


PPN = "1"
QUEUE = "D"
NODES = "1"
TIME_LIMIT = "01:00:00"
CLUSTER = "cluster-001"


# ### Initialize the endpoints

# In[ ]:


job_endpoints = JobEndpoints(*ENDPOINT_ARGS)
project_endpoints = ProjectEndpoints(*ENDPOINT_ARGS)
material_endpoints = MaterialEndpoints(*ENDPOINT_ARGS)
workflow_endpoints = WorkflowEndpoints(*ENDPOINT_ARGS)
bank_workflow_endpoints = BankWorkflowEndpoints(*ENDPOINT_ARGS)
property_endpoints = PropertiesEndpoints(*ENDPOINT_ARGS)


# Retrieve the owner and project IDs as they are needed by the endpoints. The default material is used to extract the owner ID. One can extract the owner ID from any other account's [entities](https://docs.mat3ra.com/entities-general/overview/).

# In[ ]:


owner_id = material_endpoints.list({"isDefault": True, "owner._id": ACCOUNT_ID})[0]["owner"]["_id"]
project_id = project_endpoints.list({"isDefault": True, "owner._id": ACCOUNT_ID})[0]["_id"]


# ### Create workflows
#
# Copy "ML: Train Model" and "Band Gap" bank workflows to the account's workflows. We use exabyte bank workflows which are identified by "systemName" field. The below can be adjusted to get the bank workflows by ID.

# In[ ]:


band_gap_workflow_id = copy_bank_workflow_by_system_name(bank_workflow_endpoints, "espresso-band-gap", owner_id)
ml_train_workflow_id = copy_bank_workflow_by_system_name(bank_workflow_endpoints, "exabyteml-ml-train-model", owner_id)


# ### Import materials
#
# Import materials from materials project.

# In[ ]:


train_materials = material_endpoints.import_from_materialsproject(
    MATERIALS_PROJECT_API_KEY, TRAIN_MATERIALS_PROJECT_IDS, owner_id
)
target_materials = material_endpoints.import_from_materialsproject(
    MATERIALS_PROJECT_API_KEY, TARGET_MATERIALS_PROJECT_IDS, owner_id
)


# ### Calculate Properties for "train materials"
#
# Create jobs for the "train materials".

# In[ ]:


compute = job_endpoints.get_compute(CLUSTER, PPN, NODES, QUEUE, TIME_LIMIT)
jobs = job_endpoints.create_by_ids(
    train_materials, band_gap_workflow_id, project_id, owner_id, JOB_NAME_PREFIX, compute
)


# Submit the jobs for execution.

# In[ ]:


for job in jobs:
    job_endpoints.submit(job["_id"])


# Monitor the jobs and print the status until they are all finished.

# In[ ]:


job_ids = [job["_id"] for job in jobs]
wait_for_jobs_to_finish(job_endpoints, job_ids)


# ### Build ML Train model
#
# Create ML Train job for the train materials.

# In[ ]:


name = "-".join((JOB_NAME_PREFIX, "train"))
material_ids = [m["_id"] for m in train_materials]
config = job_endpoints.get_config(material_ids, ml_train_workflow_id, project_id, owner_id, name, compute, True)
job = job_endpoints.create(config)


# Submit the train job for execution.

# In[ ]:


job_endpoints.submit(job["_id"])


# Monitor the job and print the status until it is done.

# In[ ]:


wait_for_jobs_to_finish(job_endpoints, [job["_id"]])


# ### Extract ML model as workflow
#
# The resulting trained model is extracted from the last unit (train with index 4) of the first job's subworkflow (ML: Train Model with index 0) and is further referred to as "ML predict workflow".

# In[ ]:


ml_predict_workflow = get_property_by_subworkflow_and_unit_indicies(
    property_endpoints, "workflow:ml_predict", job, 0, 4
)["data"]
ml_predict_workflow_id = ml_predict_workflow["_id"]


# Print ML predict workflow

# In[ ]:


display_JSON(ml_predict_workflow)


# ### Predict property using the model
#
# Create ML Predict job for the predict materials.

# In[ ]:


name = "-".join((JOB_NAME_PREFIX, "predict"))
material_ids = [m["_id"] for m in target_materials]
config = job_endpoints.get_config(material_ids, ml_predict_workflow_id, project_id, owner_id, name, compute, True)
job = job_endpoints.create(config)


# Submit the predict job for execution.

# In[ ]:


job_endpoints.submit(job["_id"])


# Monitor the job and print the status until its done.

# In[ ]:


wait_for_jobs_to_finish(job_endpoints, [job["_id"]])


# ### Extract predicted properties
#
# Predicted properties are extracted from the last unit (score with index 3) of the first job's subworkflow (ml_predict_subworkflow with index 0).

# In[ ]:


predicted_properties = get_property_by_subworkflow_and_unit_indicies(
    property_endpoints, "predicted_properties", job, 0, 3
)["data"]["values"]


# ### Flatten results
#
# The below for-loop iterates over the results and flatten them to form the final Pandas dataFrame.

# In[ ]:


table = []
for exabyte_id, properties in predicted_properties.items():
    material = next((m for m in target_materials if m["exabyteId"] == exabyte_id))
    band_gaps = next((v for v in properties if v["name"] == "band_gaps"))
    direct_gap = next((v for v in band_gaps["values"] if v["type"] == "direct"))["value"]
    indirect_gap = next((v for v in band_gaps["values"] if v["type"] == "indirect"))["value"]
    table.append(
        [material["_id"], material["name"], material["formula"], material["exabyteId"], direct_gap, indirect_gap]
    )


# ### Ouput results
#
# Create and print the final table as Pandas dataFrame.

# In[ ]:


headers = ["ID", "NAME", "FORMULA", "EXABYTE-ID", "DIRECT-GAP", "INDIRECT-GAP"]
df = pd.DataFrame(data=table, columns=headers)
html = dataframe_to_html(df)
html

