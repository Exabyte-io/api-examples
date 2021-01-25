#!/usr/bin/env python
# coding: utf-8

# # Get-File-From-Job
# 
# This example demonstrates how to use Exabyte RESTful API to check for and acquire files from jobs which have been run.
# This example assumes that the user is already familiar with the [creation and submission of jobs](
# create_and_submit_jobs.ipynb) using our API.
# 
# > <span style="color: orange">**IMPORTANT NOTE**</span>: In order to run this example in full, an active Exabyte.io
# account is required. Alternatively, Readers may substitute the workflow ID below with another one (an equivalent one for
# VASP, for example) and adjust extraction of the results ("Viewing job files" section). RESTful API credentials shall be
# updated in [settings](../settings.py).
# 
# 
# ## Steps
# 
# After working through this notebook, you will be able to:
# 
# 1. Import [the structure of Si](https://materialsproject.org/materials/mp-149/) from Materials Project
# 2. Set up and run a single-point calculation using Quantum Espresso.
# 3. List files currently in the job's directory
# 4. Check metadata for every file (modification date, size, etc)
# 5. Access file contents directly and print them to console
# 6. Download files to your local machine
# 
# ## Pre-requisites
# 
# The explanation below assumes that the reader is familiar with the concepts used in Exabyte platform and RESTful API. We
# outline these below and direct the reader to the original sources of information:
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
from utils import wait_for_jobs_to_finish, get_property_by_subworkow_and_unit_indicies, dataframe_to_html, \
    ensure_packages_are_installed

ensure_packages_are_installed()

# Relevant functions from the API client
from exabyte_api_client.endpoints.jobs import JobEndpoints
from exabyte_api_client.endpoints.projects import ProjectEndpoints
from exabyte_api_client.endpoints.materials import MaterialEndpoints
from exabyte_api_client.endpoints.bank_workflows import BankWorkflowEndpoints
from exabyte_api_client.endpoints.raw_properties import RawPropertiesEndpoints

# ### Create and submit the job
# 
# For this job, we'll use the workflow located [here](https://platform.exabyte.io/analytics/workflows/84DAjE9YyTFndx6z3).
# 
# This workflow is a single-point total energy calculation using Density-Functional Energy as-implemented in Quantum
# Espresso version 5.4.0.
# 
# The PBE functional is used in conjunction with an ultrasoft pseudopotential and a planewave basis set.
# 
# The material we will investigate is elemental [Silicon](https://materialsproject.org/materials/mp-149/), as-is from
# Materials Project.
# 
# > <span style="color: orange">Note</span>: This cell uses our API to copy the unit cell of silicon from Materials Project
# into your account. It then copies a workflow to get the total energy of a system using Quantum Espresso to your account.
# Finally, a job is created using the Quantum Espresso workflow for the silicon unit cell, and the job is submitted to the
# cluster. For more information, please refer to our [run-simulation-and-extract-properties](
# ./run-simulations-and-extract-properties.ipynb) notebook, located in this directory.

# In[]:


# Get some account information
project_endpoints = ProjectEndpoints(*ENDPOINT_ARGS)
project_metadata = project_endpoints.list({"isDefault": True, "owner._id": ACCOUNT_ID})[0]
project_id = project_metadata['_id']
owner_id = project_metadata['owner']['_id']

# Get a workflow for the job from the bank, and copy it to our account
bank_workflow_endpoints = BankWorkflowEndpoints(*ENDPOINT_ARGS)
BANK_WORKFLOW_ID = "84DAjE9YyTFndx6z3"
workflow_id = bank_workflow_endpoints.copy(BANK_WORKFLOW_ID, owner_id)["_id"]

# Get materials for the job
material_endpoints = MaterialEndpoints(*ENDPOINT_ARGS)
material_project_id = ["mp-149"]  # The importer expects a list
materials = material_endpoints.import_from_materialsproject(MATERIALS_PROJECT_API_KEY, material_project_id, owner_id)

# Create the job
job_endpoints = JobEndpoints(*ENDPOINT_ARGS)
job = job_endpoints.create_by_ids(materials=materials,
                                  workflow_id=workflow_id,
                                  project_id=project_id,
                                  owner_id=owner_id,
                                  prefix="Test_Job_Output")[0]

# Submit the job
job_endpoints.submit(job['_id'])
wait_for_jobs_to_finish(job_endpoints, [job['_id']])

# Monitor the jobs and print the status until they are all finished.

# ## Viewing job files
# ### Retreive a list of job files
# 
# Here, we'll get a list of all files that belong to the job.

# In[]:


files = job_endpoints.list_files(job['_id'])
paths = [file['key'] for file in files]
for path in paths:
    if 'outdir' not in path:
        print(path)

# ### Get metadata for the Output File
# The .out file is where Quantum Espresso shows its work and prints its results, so you most likely will want to view this
# files. Let's print out some of its metadata.
# 
# You'll find that we get a lot of data describing the file and its providence. Brief explanations of each entry are:
# - Key - Path to the file on the cluster
# - size - Size of the file, in bytes.
# - Bucket - The name of the cluster which ran the job.
# - Region - Which server region was used to run the job.
# - Provider - The cluster provider for the compute resources (in our case, we used AWS).
# - lastModified - Unix timestamp representing when the file was last modified.
# - name - The filename.
# - signedUrl - This is a link which can be used to download the file for a short amount of time.

# In[]:


for file in files:
    if file['name'] == 'pw_scf.out':
        output_file_metadata = file
JSON(output_file_metadata)

# ### Display file contents to console
# 
# The signedUrl gives us a place to access the file and download it. Let's read it into memory, and print out the last few
# lines of our job.

# In[]:


server_response = urllib.request.urlopen(output_file_metadata['signedUrl'])
output_file_bytes = server_response.read()

# The server returns us a bytes-string. That's useful for things like binaries or other non-human-readable data,
# but this should be decoded if we're planning to write to console.
# Because this is a human-readable text file, we'll decode it to UTF-8.
output_file = output_file_bytes.decode(encoding="UTF-8")

# Tail the last 90 lines
lines = output_file.split("\n")
for line in lines[-90:]:
    print(line)

# ### Save the input file and output file to disk.
# 
# Now that we've verified the job is done, let's go ahead and save it and its input to disk.

# In[]:


# We've already got an output file, so et's grab the input file we sent to Quantum Espresso
for file in files:
    if 'pw_scf.in' == file['name']:
        input_file_metadata = file
server_response = urllib.request.urlopen(input_file_metadata['signedUrl'])
input_file_bytes = server_response.read()

# In[]:


# Let's write the input file to disk. Note that we get files as a bytes string from the server, which is convenient for
# binaries, images, and other non-human-readable data.
# Although we could decode before writing to disk, we can just write it directly with the "wb" (write bytes) file mode.
with open(input_file_metadata['name'], 'wb') as file_descriptor:
    file_descriptor.write(input_file_bytes)

# In[]:


# Now, let's write our output file to the disk. Note that because we already decoded it, we can just use the 'w' file mode.
with open(output_file_metadata['name'], 'w') as file_descriptor:
    file_descriptor.write(output_file)
