#!/usr/bin/env python
# coding: utf-8

# In[]:


import os, sys
import string

import ase.io
import ase.neighborlist
import ase.constraints

import numpy as np

# Import settings file and utils
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)
from utils import ensure_packages_are_installed
ensure_packages_are_installed()
from material_utils import download_contcar
from settings import ENDPOINT_ARGS, ORGANIZATION_ID

# Import relevant portions of the API client
from exabyte_api_client.endpoints.jobs import JobEndpoints
from exabyte_api_client.endpoints.projects import ProjectEndpoints
from exabyte_api_client.endpoints.materials import MaterialEndpoints
from exabyte_api_client.endpoints.workflows import WorkflowEndpoints
from exabyte_api_client.endpoints.bank_workflows import BankWorkflowEndpoints


# # Query the Platform for the Slab
# 
# Get slab, acquire science

# In[]:


job_endpoint = JobEndpoints(*ENDPOINT_ARGS)
material_endpoint = MaterialEndpoints(*ENDPOINT_ARGS)

# ToDo: Refactor class and add documentation
class JobData():
    def __init__(self, jobId, job_endpoint, material_endpoint):
        self.job_endpoint = job_endpoint
        self.material_endpoint = material_endpoint
        
        self.jobId = jobId
        
        # Run the API call and get the job's name
        self.job_json = self.job_endpoint.get(jobId)
        self.name = self.job_json["name"]
        
        # Process the job name
        _, _, miller_index, end_string = self.name.split("_")
        termination_string, symbol = end_string.split(" ")
        self.termination = int("".join([char if char in string.digits else "" for char in termination_string]))
        self.miller_index = [int(index) for index in miller_index]
        self.symbol = symbol
        
        # Get the job's material
        material_id = self.job_json["_material"]["_id"]
        self.material = self.material_endpoint.get(material_id)
        
    def __repr__(self):
        return(str(self.name))
        
jobIds = (
    "tLvgnymbFi8SLNxgK",
    "5jLSJ6fJE2DPY9rqv",
    "QQbZ59P4BJnD6gRjY",
    "bvLefeo2EnKY8Lr6y",
    "wpBKynme2QGs4RqPF",
    "LZwat95GwWc4FMLFg",
    "QyYDbwjR5dhWoAvud",
    "kWfKC5mCsCB2WsovL",
    "YcuL2J4LSG56Zduqp",
    "TLzz58xDnwcyu77fX",
    "RxdryL4RJvLR5tRSF",
    "uLDz8w9FG75WGtGfp",
    "HcGbEt6ZqJvytfSac",
)
        
# Loop over list of jobIds, and generate our slabData objects
jobs = [JobData(jobId, job_endpoint, material_endpoint) for jobId in jobIds]


# # 

# # Download a Slab from the Platform
# 
# Now that we have our slab jobs, and a nice object to hold them, let's get them into PyMatGen.

# ## Get an ASE oject

# In[]:


for job in jobs:
    job.filename = f"{'_'.join(map(str, job.miller_index))}_term{job.termination}.vasp"
    download_contcar(job.jobId, job.job_endpoint, job.filename)
    job.structure = ase.io.read(job.filename)


# ## Generate the PyMatGen Object

# In[]:


# Ensure we've got at least a 5x5 angstrom surface, to help combat lateral interactions in the adsorbates
for job in jobs:
    a, b, c = job.structure.cell.lengths()
    while a < 5 or b < 5:
        if a < 5:
            job.structure *= [2,1,1]

        if b < 5:
            job.structure *= [1,2,1]
        
        a, b, c = job.structure.cell.lengths()


# ## Find the Surface Sites

# In[]:


import ase.neighborlist

def get_surface_indices(ase_structure, cns):
    surface_indices = []
    a, b, c = ase_structure.cell.lengths()
    for index, (atom, cn) in enumerate(zip(ase_structure, cns)):
        # Slabs are centered, so we can easily filter out one half of the atoms
        if cn < 12 and atom.position[2] > c/2:
            surface_indices.append(index)
    return surface_indices

for job in jobs:
    src, dest = ase.neighborlist.neighbor_list("ij", job.structure, 3)
    job.bonds = tuple(zip(src, dest))  
    job.cns = np.bincount(src)
    job.surface_indices = get_surface_indices(job.structure, job.cns)


# # CO Placement
# 
# Place CO on surface, surface norm is in Z-direction

# In[]:


def place_co_molecule(job, site):
    site_coords = job.structure[site].position
    carbon_atom = ase.Atom("C", site_coords + [0, 0, 1.43])
    oxygen_atom = ase.Atom("O", site_coords + [0, 0, 2*1.43])
    carbon_monoxide = ase.Atoms((carbon_atom, oxygen_atom))
    job.adsorption_states.append(job.structure + carbon_monoxide)
    
for job in jobs:
    job.adsorption_states = []
    for site in job.surface_indices:
        place_co_molecule(job, site)


# And finally, we'll freeze the surface to reduce the computational burden

# In[]:


def fix_all_slab_atoms(adsorption_state):
    constraint = ase.constraints.FixAtoms(indices = [atom.index for atom in adsorption_state if atom.symbol not in ("C", "O")])
    adsorption_state.set_constraint(constraint)

for job in jobs:
    for adsorption_state in job.adsorption_states:
        fix_all_slab_atoms(adsorption_state)


# # Adsorption Study
# 
# Eads=Ecomplex-(Eslab+Emolecule)

# In[]:


for job in jobs:
    for index, adsorption_state in enumerate(job.adsorption_states):
        job.adsorption_state_names = []
        adsorption_state_name = job.filename.replace(".vasp", f"site_{index}")
        job.adsorption_state_names.append(adsorption_state_name) 
        ase.io.write(adsorption_state_name + ".vasp", adsorption_state)


# ## Create the Workflow

# In[]:


# Start by finding the slab optimization workflow from the previous webinar
bank_workflow_id = "ByGrAAkeGiopxJdyu"

# Copy over the workflow
bank_workflow_endpoint = BankWorkflowEndpoints(*ENDPOINT_ARGS)
workflow = bank_workflow_endpoint.copy(bank_workflow_id, account_id=ORGANIZATION_ID)

# Copy in our constraints
vasp_unit = workflow["subworkflows"][0]["units"][0]
for input_file in vasp_unit["input"]:
    if input_file["name"] == "POSCAR":
        input_file["content"] = "{{ input.POSCAR_WITH_CONSTRAINTS }}"

# Set the names to something easy to recognize, and upload
workflow_endpoint = WorkflowEndpoints(*ENDPOINT_ARGS)
workflow['name'] = 'Constrained Adsorbate_Relaxation'
workflow = workflow_endpoint.create(workflow)


# ## Upload the Materials

# In[]:


slab_materials_set = material_endpoint.create_set({"name" : "Webinar_Cu_CO_Slabs",
                                                            "owner": {"_id": ORGANIZATION_ID}})
for job in jobs:
    job.slab_ids = []
    for jobname in job.adsorption_state_names:
        with open(jobname + ".vasp", "r") as inp:
            content = "".join(inp.readlines())
        material_json = material_endpoint.import_from_file(name=jobname, content=content, owner_id=ORGANIZATION_ID)
        job.slab_ids.append(material_json["_id"])
        material_endpoint.move_to_set(material_json["_id"], "", slab_materials_set["_id"])


# ## Submit the Jobs

# In[ ]:


projects_endpoint = ProjectEndpoints(*ENDPOINT_ARGS)
project_id = projects_endpoint.list({"isDefault": True,
                                             "owner._id": ORGANIZATION_ID})[0]["_id"]      
        
for job in jobs:
    job.job_ids = []
    for jobname, adsorption_state, material_id in zip(job.adsorption_state_names,
                                                      job.adsorption_states,
                                                      job.slab_ids):
                # First, set up the compute parameters
        # Add an extra node if there are a lot of atoms in the cell
        n_nodes = int(np.rint(len(adsorption_state)))
        # Restrict the nodes to at least 1, but no more than 2
        n_nodes = min(max(n_nodes, 1), 2)

        job_config = {"ppn": 16,
                      "queue": "OF",
                      "nodes": n_nodes,
                      "time_limit": "12:00:00",
                      "cluster": "cluster-007"}
        compute = job.job_endpoint.get_compute(**job_config)
        
        workflow_id = workflow["_id"]
        material = material_endpoint.get(material_id)
        adsorbate_job = job.job_endpoint.create_by_ids([material],
                                                       workflow_id,
                                                       project_id,
                                                       ORGANIZATION_ID,
                                                       jobname,
                                                       compute)
        job.job_ids.append(adsorbate_job[0]["_id"])

    for job_id in job.job_ids:
        job.job_endpoint.submit(job_id)


# In[ ]:




