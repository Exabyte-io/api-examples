#!/usr/bin/env python
# coding: utf-8

# # Overview
# 
# This example has been created as part of our Advanced Topics Webinar, scheduled for March 26, 2021. This webinar focused on the high-throughput screening of surfaces using Exabyte.
# 
# In this notebook, a structure for γ-Al2O3, an industrially-relevant material used in a variety of catalytic applications, both as a catalyst and as a a support. One such application is the [Claus Process](https://en.wikipedia.org/wiki/Claus_process), which is catalyzed by γ-Al2O3 and produced over 64-million tons of sulfur in 2005!

# # Execution
# 
# > <span style="color: orange">**NOTE**</span>: In order to run this example, an active Exabyte.io account is required. RESTful API credentials shall be updated in [settings](../settings.py). The generation of the credentials is also explained therein.
# 
# In addition, this notebook demonstrates the use of an "Organization ID," which allows multiple users who are part of the same organization to collaborate.
# 
# ## Import packages

# In[]:


import urllib
import os, sys

import ase.io
import ase.build.tools
import ase.cluster
from ase.visualize import view

import pymatgen.symmetry.analyzer

# Import settings file and utils
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)
from utils.generic import ensure_packages_are_installed, display_JSON, save_files
ensure_packages_are_installed()
from utils.material import get_all_slabs_and_terms, freeze_center_bulk
from settings import ENDPOINT_ARGS, ORGANIZATION_ID

# Import relevant portions of the API client
from exabyte_api_client.endpoints.jobs import JobEndpoints
from exabyte_api_client.endpoints.projects import ProjectEndpoints
from exabyte_api_client.endpoints.materials import MaterialEndpoints
from exabyte_api_client.endpoints.workflows import WorkflowEndpoints
from exabyte_api_client.endpoints.bank_workflows import BankWorkflowEndpoints


# # Visualize the Unit Cell
# 
# We begin by finding the following unit cell:
# - Gamma Al2O3, whose unit cell is in a POSCAR we constructed based on Digne, M.; Sautet, P.; Raybaud, P.; Euzen, P.; Toulhoat, H. Use of DFT to achieve a rational understanding of acid-basic properties of γ-alumina surfaces. J Catal 2004, 226, 54-68.
# 
# For this example, we have already taken the unit cell information from the literature and created a POSCAR file from it.

# In[]:


# Get gamma Al2O3 from the local disk
al2o3_poscar_filename = "../assets/gamma_alumina_digne_et_al.poscar"
al2o3_ase = ase.io.read(al2o3_poscar_filename)

# View the crystal in ASE's built-in x3d viewer
# The viewer is interactive. Press "I" to center the cell followed by "A" to zoom out
view(al2o3_ase, viewer='x3d')


# # Optimize the Unit Cells

# The next few steps will focus on optimizing the unit cell with our chosen DFT methodology, before we cleave it.
# 
# ## Create the Material
# Now that we have a unit cell, we can upload it to the platform.

# In[]:


# Create our Materials Endpoint
exabyte_materials_endpoint = MaterialEndpoints(*ENDPOINT_ARGS)

#Upload Al2O3 to the platform
with open(al2o3_poscar_filename, "r") as inp:
    content = inp.read()
    al2o3_cell_material = exabyte_materials_endpoint.import_from_file(name="Gamma_Alumina_Unit_Cell",
                                                                      content=content,
                                                                      owner_id = ORGANIZATION_ID)


# ## Create the Workflow
# Create a workflow based on the Variable-cell Relaxation workflow available [here](https://platform.exabyte.io/analytics/workflows/NAdKjws8qieKWeYnL)

# In[]:


# Start by finding a unit cell optimization workflow on the Exabyte platform.
# The bank workflow ID can be found in the bank workflows URL
exabyte_bank_workflows_endpoint = BankWorkflowEndpoints(*ENDPOINT_ARGS)
bank_workflow_id = "NAdKjws8qieKWeYnL"
al2o3_workflow = exabyte_bank_workflows_endpoint.copy(bank_workflow_id, account_id = ORGANIZATION_ID)

# Set the names to something easy to recognize, and upload
exabyte_workflows_endpoint = WorkflowEndpoints(*ENDPOINT_ARGS)
al2o3_workflow['name'] = "Alumina_Cell_Relax"
al2o3_cellopt_workflow = exabyte_workflows_endpoint.create(al2o3_workflow)


# ## Create the Job
# Then create the jobs to optimize these unit cells

# In[]:


# Get the default projectID for the organization account (this can be substituted with a userID)
exabyte_projects_endpoint = ProjectEndpoints(*ENDPOINT_ARGS)
project_id = exabyte_projects_endpoint.list({"isDefault": True,
                                             "owner._id": ORGANIZATION_ID})[0]["_id"]

# Create the compute configuration for the jobs on Azure
exabyte_jobs_endpoint = JobEndpoints(*ENDPOINT_ARGS)
job_config = {"ppn": 16,
              "queue": "OF",
              "nodes": 1,
              "time_limit": "00:20:00",
              "cluster": "cluster-007"}
compute = exabyte_jobs_endpoint.get_compute(**job_config)

# Create the Al2O3 job
al2o3_job = exabyte_jobs_endpoint.create_by_ids([al2o3_cell_material],
                                                al2o3_cellopt_workflow["_id"],
                                                project_id,
                                                ORGANIZATION_ID,
                                                "al2o3_cellopt",
                                                compute)


# Finally, submit the jobs

# In[]:


exabyte_jobs_endpoint.submit(al2o3_job[0]["_id"])


# # Generate the Slabs

# We begin by extracting the relaxed unit cell from the jobs

# In[]:


# Extract job ID
al2o3_job_id = al2o3_job[0]["_id"]

save_files(al2o3_job_id, exabyte_jobs_endpoint, "CONTCAR", "al2o3_relaxed.vasp")


# The below code uses Pymatgen to find all unique planes in a crystal. It then generates slabs with every possible termination in that plane. Finally, asymmetric slabs are filtered out, and the slabs are saved to a dictionary. The dictionary's format is: `{miller-index: {termination: {"slab": slab} }`. Note that here, to take advantage of several convenient functions in ASE, we have the function output ASE Atoms objects instead of PyMatGen objects.

# In[]:


al2o3_relaxed_pymatgen = pymatgen.core.structure.Structure.from_file("al2o3_relaxed.vasp")

al2o3_slabs = get_all_slabs_and_terms(al2o3_relaxed_pymatgen, thickness=3, is_by_layers=True)


# # Write the Slabs to Disk
# 
# To finish preparing our structures for VASP, we do the following for both sets of slabs:
# 1. Center the slab, and adjust the vacuum to 10Å
# 2. Freeze the center layer of the slab to simulate the bulk
# 3. Write the slab to a vasp-format file

# In[]:


# Declare the vacuum size we want to have (for easy adjustment by future users)
vacuum_size = 10
for miller_index, term_dict in al2o3_slabs.items():
    for term, surface in term_dict.items():
        # Sort the slab by symbol, to help make the POTCAR file smaller
        slab = ase.build.tools.sort(surface["slab"])

        # Center the slab's coordinates and adjust the vacuum to vacuum_size
        # Note that the "vacuum" argument refers to the amount of vacuum on either side
        # So, we divide by 2
        slab.center(vacuum=vacuum_size/2, axis=2)

        # Add a FixAtoms constraint to freeze the center layer of atoms
        freeze_center_bulk(slab)

        # Write the slab to a file
        formula = slab.get_chemical_formula(empirical=True)
        filename = f"{formula}_{miller_index}_term{term}.vasp"
        ase.io.write(filename, slab)


# # Upload POSCARS
# 
# Now that we have a set of prepared POSCARs, we can upload them to our user account on Exabyte. Let's create a couple material sets to hold these.

# In[]:


al2o3_materials_set = exabyte_materials_endpoint.create_set({"name" : "Webinar_Slabs_Al2O3",
                                                             "owner" : {"_id" : ORGANIZATION_ID}})


# And finally, we will upload all of the slabs to the platform

# In[ ]:


for miller_index, term_dict in al2o3_slabs.items():
    for term, surface in term_dict.items():
        slab = surface["slab"]

        # Get the filename we generated above
        formula = slab.get_chemical_formula(empirical=True)
        material_name = f"{formula}_{miller_index}_term{term}"
        filename = f"{material_name}.vasp"

        # Import the material to Exabyte, and place it in the correct material set
        with open(filename, "r") as inp:
            content = inp.read()
        material_json = exabyte_materials_endpoint.import_from_file(name=material_name, content=content,
                                                                    owner_id=ORGANIZATION_ID)
        exabyte_materials_endpoint.move_to_set(material_json["_id"], "", al2o3_materials_set["_id"])


# In[ ]:




