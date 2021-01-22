#!/usr/bin/env python
# coding: utf-8

# # Get-File-From-Job
# 
# This example demonstrates how to use Exabyte RESTful API to check for and acquire files from jobs which have been run. This example assumes that the user is already familiar with the [creation and submission of jobs](create_and_submit_jobs.ipynb) using our API.
# 
# > <span style="color: orange">**IMPORTANT NOTE**</span>: In order to run this example in full, an active Exabyte.io account with access to VASP (Vienna ab-initio simulations package) is required. Alternatively, Readers may substitute the workflow ID below with another one (an equivalent one for Quantum ESPRESSO, for example) and adjust extraction of the results ("Extract results" section). RESTful API credentials shall be updated in [settings](../settings.py).
# 
# 
# ## Steps
# 
# In this notebook, we will accomplish the following:
# 
# 1. Import the structure of cubic Ge from Materials Project
# 2. Set up and run a single-point calculation using VASP
# 3. List the files currently in the job's directory
# 4. Get download links for the OUTCAR and CONTCAR files of a VASP run
# 5. Decode the OUTCAR file into UTF-8, and write the last few contents to the screen
# 6. Save both the OUTCAR and CONTCAR files to disk.
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


from IPython.display import JSON, display
import os
import sys
import urllib

# Import settings file and utils file
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path: sys.path.append(module_path)
from settings import ENDPOINT_ARGS, ACCOUNT_ID, MATERIALS_PROJECT_API_KEY
from utils import wait_for_jobs_to_finish, get_property_by_subworkow_and_unit_indicies, dataframe_to_html, ensure_packages_are_installed
ensure_packages_are_installed()

# Relevant functions from the API client
from exabyte_api_client.endpoints.jobs import JobEndpoints
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


MATERIALS_PROJECT_ID = ["mp-32"]
MATERIALS_SET_NAME = "materials-set"
TAGS = ["tag1", "tag2"]


# #### Jobs
# 
# Parameters for the jobs to be ran for the imported materials:
# 
# - **JOB_NAME_PREFIX**: prefix to be used for the job name with "{JOB_NAME_PREFIX} {FORMULA}" convention (e.g.  "Job Name Prefix - SiGe")
# - **JOBS_SET_NAME**: the name of the jobs set

# In[]:


JOB_NAME_PREFIX = "Ge_Test"
JOBS_SET_NAME = "jobs-set"


# #### Workflow
# 
# This example is based on [this](https://platform.exabyte.io/analytics/workflows/SnGspffD7nPsBHJaQ) bank workflow which is later copied to the account workflows collection.

# In[]:


BANK_WORKFLOW_ID = "SnGspffD7nPsBHJaQ"


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


materials = material_endpoints.import_from_materialsproject(MATERIALS_PROJECT_API_KEY, MATERIALS_PROJECT_ID, owner_id)


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


# ### Submit jobs
# Submit the jobs for execution.

# In[]:


for job in jobs: job_endpoints.submit(job["_id"])


# Monitor the jobs and print the status until they are all finished.

# In[]:


job_ids = [job["_id"] for job in jobs]
wait_for_jobs_to_finish(job_endpoints, job_ids)


# ### Retreive a list of job files
# 
# Here, we'll get a list of all files that belong to the job.

# In[]:


job = jobs[0]
files = job_endpoints.list_files(job['_id'])
filenames = [file['key'] for file in files]
for filename in filenames:
    print(filename)


# # Get more information on the OUTCAR
# The OUTCAR file is where VASP shows it work, and contains a great deal of useful information that you might want to parse programatically. Very commonly in computational chemistry, users will "tail" the OUTCAR (that is, only printing the last several lines) to check whether it has completed.
# 
# Focusing on the OUTCAR, let's see what metadata is stored alongside a file.

# In[]:


for file in files:
    if file['name'] == 'OUTCAR':
        outcar_metadata = file
    elif file['name'] == 'CONTCAR':
        contcar_metadata = file
        
JSON(outcar_metadata)


# ### Check the OUTCAR
# 
# We get a lot of data describing the file and its providence. Brief explanations of each entry are:
# - Key - Path to the file on the cluster
# - size - Size of the file, in bytes.
# - Bucket - The name of the cluster which ran the job.
# - Region - Which server region was used to run the job.
# - Provider - The cluster provider for the compute resources (in our case, we used AWS).
# - lastModified - Unix timestamp representing when the file was last modified.
# - name - The filename.
# - signedUrl - This is a link which can be used to download the file for a short amount of time.
# 
# Below, we'll read data into memory from the signedUrl. This could be used to write the file to the disk, but for now we'll just "tail" the last few lines of the file.

# In[]:


server_response = urllib.request.urlopen(outcar_metadata['signedUrl'])
outcar = server_response.read()

# The outcar is currently stored as a bytes string in Python. That's useful for things like binaries or other non-human-readable data, but this should be decoded if it's intended to be human-readable.
# Because this is a human-readable text file, we'll decode it to UTF-8.
outcar = outcar.decode(encoding="UTF-8")

# Tail the last 28 lines
lines = outcar.split("\n")
for line in lines[-28:]:
    print(line)


# ### Write the OUTCAR and CONTCAR to disk
# 
# Now that the VASP job is done, we'll write the OUTCAR (contaiing a detailed history of the calculation) along with the CONTCAR (which, in the case of optimization jobs, would contain a converged structure).

# In[]:


server_response = urllib.request.urlopen(contcar_metadata['signedUrl'])
contcar = server_response.read()

# Recall that earlier, we decoded the OUTCAR to UTF-8 format. So, "w" as the file speification will be fine here.
with open("OUTCAR", "w") as outp:
    outp.write(outcar)

# To demonstrate that files don't require decoding if they're just being downloaded, we'll use the "wb" file specification (short for "write bytes") to save the CONTCAR.
with open("CONTCAR", "wb") as outp:
    outp.write(contcar)
