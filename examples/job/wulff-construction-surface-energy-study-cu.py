#!/usr/bin/env python
# coding: utf-8

# # Overview
# 
# This example has been created as part of our Advanced Topics Webinar, scheduled for March 26, 2021. This webinar focused on the high-throughput screening of surfaces using Exabyte.
# 
# In this notebook, we investigate the stability of a variety of different surfaces of Cu. After the surface energies are all automatically calculated, we finish by creating a Wulff Construction of Cu.

# In[]:


import os, sys

import ase.io
import ase.cluster

import pymatgen.ext.matproj
import pymatgen.io.ase
import pymatgen.symmetry.analyzer
import numpy as np

# Import settings file and utils
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)
from utils.generic import ensure_packages_are_installed, save_files
ensure_packages_are_installed()
from utils.material import get_all_slabs_and_terms
from utils.material import get_vasp_total_energy, get_slab_area, get_surface_energy
from settings import MATERIALS_PROJECT_API_KEY, ENDPOINT_ARGS, ORGANIZATION_ID

# Import relevant portions of the API client
from exabyte_api_client.endpoints.jobs import JobEndpoints
from exabyte_api_client.endpoints.projects import ProjectEndpoints
from exabyte_api_client.endpoints.materials import MaterialEndpoints
from exabyte_api_client.endpoints.workflows import WorkflowEndpoints
from exabyte_api_client.endpoints.bank_workflows import BankWorkflowEndpoints


# # Get Unit Cells
# 
# We begin by finding the following unit cell:
# - Cu, whose conventional unit cell is taken directly from Materials Project

# In[]:


def get_best_structure_by_formula(formula: str) -> pymatgen.ext.matproj.ComputedStructureEntry:
    # Create the Materials Project API object
    materials_project_api = pymatgen.ext.matproj.MPRester(MATERIALS_PROJECT_API_KEY)

    # Query for possible cells of Cu
    cu_structs = materials_project_api.get_entries(formula,
                                                   conventional_unit_cell=True,
                                                   sort_by_e_above_hull=True,
                                                   inc_structure="initial")
    best_cu_struct = cu_structs[0]
    return best_cu_struct

# Get Cu from Materials Project
cu_pymatgen = get_best_structure_by_formula("Cu").structure
cu_ase = pymatgen.io.ase.AseAtomsAdaptor.get_atoms(cu_pymatgen)


# Next, upload the unit cell to the Exabyte.IO Platform.
# 
# Note that this is in contrast to the way we import directly from Materials Project via IDs as in `materials/import_from_materials_project.ipynb` or in `material/api_interoperability_showcase.ipynb`. This is done deliberately to demonstrate the flexibility offered by the API.

# In[]:


# Create our Materials Endpoint
exabyte_materials_endpoint = MaterialEndpoints(*ENDPOINT_ARGS)

# Uplad Cu
cu_poscar_filename = "cu.vasp"
ase.io.write(cu_poscar_filename, cu_ase)
with open(cu_poscar_filename, "r") as inp:
    content = inp.read()
    cu_cell_material = exabyte_materials_endpoint.import_from_file(name="Copper_Unit_Cell", 
                                                                   content=content,
                                                                   owner_id=ORGANIZATION_ID)


# # Optimize the Unit Cell
# 
# We begin by specifying our INCAR file

# In[]:


# INCAR based on https://www.nature.com/articles/sdata201680#Sec2
incar_content = ["EDIFF = 1e-6",
                 "EDIFFG = -0.01",
                 "ISMEAR = 2",
                 "IBRION = 2",
                 "ISIF = 3",
                 "KPAR = 4",
                 "NSW = 300",
                 "ENCUT = 400"]

kpoints_content = ["Automatic mesh",
                   "0",
                   "Gamma",
                   "8 8 8",
                   "0 0 0"]


# Create a workflow based on the Variable-cell Relaxation workflow available [here](https://platform.exabyte.io/analytics/workflows/NAdKjws8qieKWeYnL)

# In[]:


# Start by finding a unit cell optimization workflow on the Exabyte platform. The bank workflow ID can be found in the bank workflows URL
bank_workflow_id = "NAdKjws8qieKWeYnL"

# Set up our Bank Workflow endpoint and copy over the workflow
exabyte_bank_workflows_endpoint = BankWorkflowEndpoints(*ENDPOINT_ARGS)
cu_workflow = exabyte_bank_workflows_endpoint.copy(bank_workflow_id,
                                                   account_id = ORGANIZATION_ID)

# Update the workflow with our custom INCAR file
vasp_unit = cu_workflow["subworkflows"][0]["units"][0]
for input_file in vasp_unit["input"]:
    if input_file["name"] == "INCAR":
        input_file['content'] = "\n".join(incar_content)
    elif input_file["name"] == "KPOINTS":
        input_file['content'] = "\n".join(kpoints_content)

# Set the names to something easy to recognize, and upload
cu_workflow['name'] = 'Copper_Cell_Relax'
exabyte_workflows_endpoint = WorkflowEndpoints(*ENDPOINT_ARGS)
cu_cellopt_workflow = exabyte_workflows_endpoint.create(cu_workflow)


# Create the job for the unit cell optimization

# In[]:


# Get the default projectID for the organization account (this can be substituted with a userID)
exabyte_projects_endpoint = ProjectEndpoints(*ENDPOINT_ARGS)
project_id = exabyte_projects_endpoint.list({"isDefault": True,
                                             "owner._id": ORGANIZATION_ID})[0]["_id"]

# Create the compute configuration for the jobs on Azure
exabyte_jobs_endpoint = JobEndpoints(*ENDPOINT_ARGS)
job_config = {"ppn": 4,
              "queue": "D",
              "nodes": 1,
              "time_limit": "00:10:00",
              "cluster": "cluster-001"}
compute = exabyte_jobs_endpoint.get_compute(**job_config)

# Create the Cu job
cu_job = exabyte_jobs_endpoint.create_by_ids([cu_cell_material],
                                             cu_cellopt_workflow["_id"],
                                             project_id,
                                             ORGANIZATION_ID,
                                             "cu_cellopt",
                                             compute)


# Finally, submit the job

# In[]:


exabyte_jobs_endpoint.submit(cu_job[0]["_id"])


# # Generate the Slabs

# We begin by extracting the relaxed unit cell from the job
# 

# In[]:


# Get a list of files for each
cu_job_id = cu_job[0]["_id"]

save_files(cu_job_id, exabyte_jobs_endpoint, "CONTCAR", "cu_relaxed.vasp")


# The below code uses Pymatgen to find all unique planes in a crystal. It then generates slabs with every possible termination in that plane. Finally, asymmetric slabs are filtered out, and the slabs are saved to a dictionary. The dictionary's format is: `{miller-index: {termination: {"slab": slab} }`. Note that here, to take advantage of several convenient functions in ASE, we have the function output ASE Atoms objects instead of PyMatGen objects.

# In[]:


cu_relaxed_pymatgen = pymatgen.core.structure.Structure.from_file("cu_relaxed.vasp")

# It is important to standardize the cell before cleaving planes
standardizer = pymatgen.symmetry.analyzer.SpacegroupAnalyzer(cu_relaxed_pymatgen)
cu_relaxed_pymatgen = standardizer.get_conventional_standard_structure()

cu_slabs = get_all_slabs_and_terms(cu_relaxed_pymatgen, thickness = 15, is_by_layers=False)


# # Write Slabs to POSCARs
# 
# A critical part of any VASP calculation is the actual POSCAR file - the file that contains the atomic coordinates.
# In this section, we take the slabs we generated above, and prepare them for the VASP calculation.

# In[]:


# Declare the vacuum size we want to have (for easy adjustment by future users)
vacuum_size = 15
for miller_index, term_dict in cu_slabs.items():
    for term, surface in term_dict.items():
        slab = surface["slab"]
        # Center the slab's coordinates and adjust the vacuum to vacuum_size
        # Note that the "vacuum" argument refers to the amount of vacuum on either side - so we divide by 2
        slab.center(vacuum=vacuum_size/2, axis=2)

        # Write the slab to a file
        formula = slab.get_chemical_formula(empirical=True)
        filename = f"{formula}_{miller_index}_term{term}.vasp"
        ase.io.write(filename, slab)


# # Upload POSCARS to Exabyte.io
# 
# Now that we have a set of prepared POSCARs, we can upload them to our user account on Exabyte. Let's create a couple material sets to hold these.

# In[]:


cu_materials_set = exabyte_materials_endpoint.create_set({"name" : "Webinar_Cu_Wulff_Slabs",
                                                         "owner": {"_id": ORGANIZATION_ID}})


# We want to keep track of material IDs as we go along, so we'll add this key to our Surface dictionary

# In[]:


for miller_index, term_dict in cu_slabs.items():
    for term, surface in term_dict.items():
        slab = surface["slab"]

        # Get the filename we generated above
        formula = slab.get_chemical_formula(empirical=True)
        material_name = f"{formula}_{miller_index}_term{term}"
        filename = f"{material_name}.vasp"

        # Import the material to Exabyte, and place it in the correct material set
        with open(filename, "r") as inp:
            content = inp.read()
        material_json = exabyte_materials_endpoint.import_from_file(name=material_name,
                                                                    content=content,
                                                                    owner_id=ORGANIZATION_ID)
        exabyte_materials_endpoint.move_to_set(material_json["_id"], "", cu_materials_set["_id"])

        # Adjust our dictionary to keep track of the Material ID
        surface["material_id"] = material_json["_id"]


# # Optimize the Cu Slabs
# 
# Now it is time to optimize the Cu slabs. We begin by specifying our INCAR file. It is largely the same as the previous INCAR, except with ISIF set to 2, to avoid changing the cell dimensions.

# In[]:


incar_content = ["EDIFF = 1e-6",
                 "EDIFFG = -0.01",
                 "ISMEAR = 2",
                 "IBRION = 2",
                 "ISIF = 2",
                 "KPAR = 4",
                 "NSW = 300",
                 "ENCUT = 400"]
kpoints_content = ["Automatic mesh",
                   "0",
                   "Gamma",
                   "8 8 1",
                   "0 0 0"]


# Next, we'll create the surface relaxation jobs. We start by cloning in a fixed cell relaxation workflow and adjusting it to suit our purposes.
# The workflow we will choose can be found [here](https://platform.exabyte.io/analytics/workflows/rTEtXntXo3ScGhi7q)

# In[]:


# Start by finding a fixed cell optimization workflow on the Exabyte platform.
bank_workflow_id = "rTEtXntXo3ScGhi7q"

# Copy over the workflow
workflow = exabyte_bank_workflows_endpoint.copy(bank_workflow_id,
                                                account_id=ORGANIZATION_ID)

# Add in our INCAR
vasp_unit = workflow["subworkflows"][0]["units"][0]
for input_file in vasp_unit["input"]:
    if input_file["name"] == "INCAR":
        input_file["content"] = "\n".join(incar_content)
    elif input_file["name"] == "KPOINTS":
        input_file["content"] = "\n".join(kpoints_content)

# Set the names to something easy to recognize, and upload
workflow['name'] = 'Copper_Slab_Relax'
workflow = exabyte_workflows_endpoint.create(workflow)


# Next, create the slab relaxation jobs

# In[]:


for miller_index, term_dict in cu_slabs.items():
    for term, surface in term_dict.items():
        # First, set up the compute parameters
        # Add an extra node if there are a lot of atoms in the cell
        n_nodes = int(np.rint(len(surface["slab"])/16))
        # Restrict the nodes to at least 1, but no more than 2
        n_nodes = min(max(n_nodes, 1), 2)

        job_config = {"ppn": 16,
                      "queue": "OF",
                      "nodes": n_nodes,
                      "time_limit": "12:00:00",
                      "cluster": "cluster-007"}
        compute = exabyte_jobs_endpoint.get_compute(**job_config)

        # Determine what we should name these jobs
        job_prefix = f"cu_slab_{miller_index}_Term{term}"

        # Get the material and workflow
        material_id = surface["material_id"]
        workflow_id = workflow["_id"]
        material = exabyte_materials_endpoint.get(material_id)

        # Set up the job and record the JOB ID, then submit the job
        cu_slab_job = exabyte_jobs_endpoint.create_by_ids([material],
                                                          workflow_id,
                                                          project_id,
                                                          ORGANIZATION_ID,
                                                          job_prefix,
                                                          compute)
        surface["job_id"] = cu_slab_job[0]["_id"]
        exabyte_jobs_endpoint.submit(surface["job_id"])


# # Extract Bulk and Slab Energy
# 
# 
# After our jobs complete, we can calculate the surface energy as follows:
# 
# (E_Slab - E_bulk * (N_Slab / N_Bulk)) / (2A)
# 
# We begin by extracting the total energy from the Cu unit cell

# In[]:


cu_job_id = cu_job[0]["_id"]
cu_bulk_energy = get_vasp_total_energy(cu_job_id, exabyte_jobs_endpoint)


# Then, we extract the total energy from each of the slab runs

# In[]:


for miller_index, term_dict in cu_slabs.items():
    for term, surface in term_dict.items():
        surface["slab_energy"] = get_vasp_total_energy(surface["job_id"], exabyte_jobs_endpoint)


# # Calculate the Surface Energy
# We'll iterate over each slab and calculate the surface energy

# In[]:


for miller_index, term_dict in cu_slabs.items():
    for term, surface in term_dict.items():
        slab = surface["slab"]
        # Slab and bulk energy
        e_slab = surface['slab_energy']
        e_bulk = cu_bulk_energy

        # N bulk and N Slab
        cu_cell = ase.io.read("cu_relaxed.vasp")
        n_bulk = len(cu_cell)
        n_slab = len(slab)

        # Slab Area
        vec_a = slab.cell[0]
        vec_b = slab.cell[1]
        area = get_slab_area(vec_a, vec_b)

        print(miller_index, e_slab, e_bulk, n_slab, n_bulk, area)

        # Surface Energy
        surface_energy = get_surface_energy(e_slab=e_slab, e_bulk=e_bulk,
                                            n_slab=n_slab, n_bulk=n_bulk,
                                            a=area)
        surface["surface_energy"] = surface_energy


# # Wulff Construction
# 
# Finally, now that we have a collection of surfaces and their energies, we can perform the Wulff Construction

# In[]:


surfaces = []
energies = []
for miller_index, term_dict in cu_slabs.items():
    # Turn the string miller index into a list, for ASE
    i,j,k = map(int, miller_index)
    plane = (i,j,k)
    surfaces.append(plane)
    # Find the best termination for the surface
    best_term_energy = np.inf
    for term, surface in term_dict.items():
        best_term_energy = min(best_term_energy, surface["surface_energy"])
    energies.append(best_term_energy)

cluster_size = 147
wulff = ase.cluster.wulff_construction("Cu", surfaces=surfaces,
                                       energies=energies, size=cluster_size,
                                       structure = "fcc")

ase.io.write("Cu_Wulff.xyz", wulff)


# In[]:


for surf, en in sorted(zip(surfaces, energies), key = lambda i: i[1]):
    print(surf, np.round(en,3))


# In[ ]:




