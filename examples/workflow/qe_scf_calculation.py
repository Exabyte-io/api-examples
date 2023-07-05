#!/usr/bin/env python
# coding: utf-8

# <a href="https://colab.research.google.com/github/Exabyte-io/api-examples/blob/dev/examples/workflow/qe_scf_calculation.ipynb" target="_parent">
# <img alt="Open in Google Colab" src="https://user-images.githubusercontent.com/20477508/128780728-491fea90-9b23-495f-a091-11681150db37.jpeg" width="150" border="0">
# </a>

# # Quantum Espresso SCF calculation via API
# 

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


# In[ ]:


from utils.settings import ENDPOINT_ARGS, ACCOUNT_ID
from utils.generic import display_JSON, wait_for_jobs_to_finish

from exabyte_api_client.endpoints.workflows import WorkflowEndpoints
from exabyte_api_client.endpoints.materials import MaterialEndpoints
from exabyte_api_client.endpoints.jobs import JobEndpoints


# In[ ]:


# Initialize a helper class to interact with WorkflowEndpoints
workflow_endpoints = WorkflowEndpoints(*ENDPOINT_ARGS)
material_endpoints = MaterialEndpoints(*ENDPOINT_ARGS)
job_endpoints = JobEndpoints(*ENDPOINT_ARGS)


# #### Create a Quantum Espresso workflow for SCF calculation

# In[ ]:


# payload for workflow creation
WORKFLOW_BODY = {
    "name": "Silicon SCF",
    "subworkflows": [
        {
            "name": "Total energy",
            "application": {
                "_id": "u2qtBhseFfQMWRYND",
                "name": "espresso",
                "shortName": "qe",
                "summary": "Quantum Espresso",
                "build": "Intel",
                "version": "6.3",
                "isDefault": False,
                "hasAdvancedComputeOptions": True,
                "schemaVersion": "2022.8.16",
                "createdAt": "2023-04-13T09:42:52.622Z",
                "createdBy": "0",
                "updatedAt": "2023-06-15T20:01:39.744Z",
                "updatedBy": "0",
            },
            "properties": [
                "total_energy",
                "total_energy_contributions",
                "pressure",
                "fermi_energy",
                "atomic_forces",
                "total_force",
                "stress_tensor",
            ],
            "model": {
                "type": "dft",
                "subtype": "gga",
                "method": {"type": "pseudopotential", "subtype": "us", "data": {"searchText": ""}},
                "functional": {"slug": "pbe"},
                "refiners": [],
                "modifiers": [],
            },
            "units": [
                {
                    "type": "execution",
                    "application": {
                        "_id": "u2qtBhseFfQMWRYND",
                        "name": "espresso",
                        "shortName": "qe",
                        "summary": "Quantum Espresso",
                        "build": "Intel",
                        "version": "6.3",
                        "isDefault": False,
                        "hasAdvancedComputeOptions": True,
                        "schemaVersion": "2022.8.16",
                        "createdAt": "2023-04-13T09:42:52.622Z",
                        "createdBy": "0",
                        "updatedAt": "2023-06-15T20:01:39.744Z",
                        "updatedBy": "0",
                    },
                    "flowchartId": "066bb68919ae51b993ce7739",
                    "status": "idle",
                    "statusTrack": [],
                    "results": [
                        {"name": "total_energy"},
                        {"name": "total_energy_contributions"},
                        {"name": "pressure"},
                        {"name": "fermi_energy"},
                        {"name": "atomic_forces"},
                        {"name": "total_force"},
                        {"name": "stress_tensor"},
                    ],
                    "monitors": [{"name": "standard_output"}, {"name": "convergence_electronic"}],
                    "preProcessors": [],
                    "postProcessors": [],
                    "head": True,
                    "input": [
                        {
                            "content": "{% if subworkflowContext.MATERIAL_INDEX %}\n{%- set input = input.perMaterial[subworkflowContext.MATERIAL_INDEX] -%}\n{% endif -%}\n&CONTROL\n    calculation = 'scf'\n    title = ''\n    verbosity = 'low'\n    restart_mode = '{{ input.RESTART_MODE }}'\n    wf_collect = .true.\n    tstress = .true.\n    tprnfor = .true.\n    outdir = {% raw %}'{{ JOB_WORK_DIR }}/outdir'{% endraw %}\n    wfcdir = {% raw %}'{{ JOB_WORK_DIR }}/outdir'{% endraw %}\n    prefix = '__prefix__'\n    pseudo_dir = {% raw %}'{{ JOB_WORK_DIR }}/pseudo'{% endraw %}\n/\n&SYSTEM\n    ibrav = {{ input.IBRAV }}\n    nat = {{ input.NAT }}\n    ntyp = {{ input.NTYP }}\n    ecutwfc = {{ cutoffs.wavefunction }}\n    ecutrho = {{ cutoffs.density }}\n    occupations = 'smearing'\n    degauss = 0.005\n/\n&ELECTRONS\n    diagonalization = 'david'\n    diago_david_ndim = 4\n    diago_full_acc = .true.\n    mixing_beta = 0.3\n    startingwfc = 'atomic+random'\n/\n&IONS\n/\n&CELL\n/\nATOMIC_SPECIES\n{{ input.ATOMIC_SPECIES }}\nATOMIC_POSITIONS crystal\n{{ input.ATOMIC_POSITIONS }}\nCELL_PARAMETERS angstrom\n{{ input.CELL_PARAMETERS }}\nK_POINTS automatic\n{% for d in kgrid.dimensions %}{{d}} {% endfor %}{% for s in kgrid.shifts %}{{s}} {% endfor %}\n",
                            "name": "pw_scf.in",
                            "contextProviders": [
                                {"name": "KGridFormDataManager"},
                                {"name": "QEPWXInputDataManager"},
                                {"name": "PlanewaveCutoffDataManager"},
                            ],
                            "applicationName": "espresso",
                            "executableName": "pw.x",
                            "rendered": "&CONTROL\n    calculation = 'scf'\n    title = ''\n    verbosity = 'low'\n    restart_mode = 'from_scratch'\n    wf_collect = .true.\n    tstress = .true.\n    tprnfor = .true.\n    outdir = '{{ JOB_WORK_DIR }}/outdir'\n    wfcdir = '{{ JOB_WORK_DIR }}/outdir'\n    prefix = '__prefix__'\n    pseudo_dir = '{{ JOB_WORK_DIR }}/pseudo'\n/\n&SYSTEM\n    ibrav = 0\n    nat = 2\n    ntyp = 1\n    ecutwfc = 40\n    ecutrho = 200\n    occupations = 'smearing'\n    degauss = 0.005\n/\n&ELECTRONS\n    diagonalization = 'david'\n    diago_david_ndim = 4\n    diago_full_acc = .true.\n    mixing_beta = 0.3\n    startingwfc = 'atomic+random'\n/\n&IONS\n/\n&CELL\n/\nATOMIC_SPECIES\nSi 28.0855 si_pbe_gbrv_1.0.upf\nATOMIC_POSITIONS crystal\nSi  0.000000000 0.000000000 0.000000000 \nSi  0.250000000 0.250000000 0.250000000 \nCELL_PARAMETERS angstrom\n3.348920236 0.000000000 1.933500000\n1.116307420 3.157392040 1.933500000\n0.000000000 0.000000000 3.867000000\nK_POINTS automatic\n5 5 5 0 0 0 \n",
                        }
                    ],
                    "context": {
                        "kgridExtraData": {"materialHash": "a665723ef7429caef6ca89385fe25bae"},
                        "kgrid": {
                            "dimensions": [5, 5, 5],
                            "shifts": [0, 0, 0],
                            "reciprocalVectorRatios": [1, 1, 1],
                            "gridMetricType": "KPPRA",
                            "gridMetricValue": 250,
                            "preferGridMetric": False,
                        },
                        "isKgridEdited": True,
                        "subworkflowContext": {},
                    },
                    "executable": {
                        "isDefault": True,
                        "hasAdvancedComputeOptions": True,
                        "postProcessors": ["remove_non_zero_weight_kpoints"],
                        "monitors": ["standard_output", "convergence_ionic", "convergence_electronic"],
                        "results": [
                            "atomic_forces",
                            "band_gaps",
                            "density_of_states",
                            "fermi_energy",
                            "pressure",
                            "stress_tensor",
                            "total_energy",
                            "total_energy_contributions",
                            "total_force",
                            "final_structure",
                            "magnetic_moments",
                            "reaction_energy_barrier",
                            "reaction_energy_profile",
                            "potential_profile",
                            "charge_density_profile",
                        ],
                        "name": "pw.x",
                    },
                    "flavor": {
                        "isDefault": True,
                        "input": [{"name": "pw_scf.in"}],
                        "results": [
                            "total_energy",
                            "total_energy_contributions",
                            "pressure",
                            "fermi_energy",
                            "atomic_forces",
                            "total_force",
                            "stress_tensor",
                        ],
                        "monitors": ["standard_output", "convergence_electronic"],
                        "applicationName": "espresso",
                        "executableName": "pw.x",
                        "name": "pw_scf",
                        "executable": {
                            "isDefault": True,
                            "hasAdvancedComputeOptions": True,
                            "postProcessors": ["remove_non_zero_weight_kpoints"],
                            "monitors": ["standard_output", "convergence_ionic", "convergence_electronic"],
                            "results": [
                                "atomic_forces",
                                "band_gaps",
                                "density_of_states",
                                "fermi_energy",
                                "pressure",
                                "stress_tensor",
                                "total_energy",
                                "total_energy_contributions",
                                "total_force",
                                "final_structure",
                                "magnetic_moments",
                                "reaction_energy_barrier",
                                "reaction_energy_profile",
                                "potential_profile",
                                "charge_density_profile",
                            ],
                            "name": "pw.x",
                        },
                    },
                    "name": "pw_scf",
                }
            ],
        }
    ],
    "units": [{"name": "Total energy", "type": "subworkflow", "head": True, "schemaVersion": "2022.8.16"}],
    "properties": [
        "total_energy",
        "total_energy_contributions",
        "pressure",
        "fermi_energy",
        "atomic_forces",
        "total_force",
        "stress_tensor",
    ],
}


# In[ ]:


# create workflow
WORKFLOW_RESP = workflow_endpoints.create(WORKFLOW_BODY)


# #### Get default material from the user account

# In[ ]:


default_material = material_endpoints.list({"isDefault": True, "owner._id": ACCOUNT_ID})[0]
material_id = default_material["_id"]


# In[ ]:


# job creation payload
JOB_BODY = {
    "name": "SCF Calculation",
    "compute": {
        "ppn": 2,
        "nodes": 1,
        "queue": "OR",
        "cluster": {"fqdn": "master-production-20160630-cluster-001.exabyte.io"},
    },
    "_project": {"slug": "pranab-default"},
    "workflow": WORKFLOW_RESP,
    "_material": {"_id": material_id},
}


# #### Create and submit job

# In[ ]:


# create job
JOB_RESP = job_endpoints.create(JOB_BODY)


# In[ ]:


# submit job
job_endpoints.submit(JOB_RESP["_id"])


# In[ ]:


wait_for_jobs_to_finish(job_endpoints, [JOB_RESP["_id"]])


# #### Get output file, perform post processing, and make plots

# In[ ]:


files = job_endpoints.list_files(JOB_RESP["_id"])
for file in files:
    if file["name"] == "pw_scf.out":
        output_file_metadata = file

import urllib

server_response = urllib.request.urlopen(output_file_metadata["signedUrl"])
output_file_bytes = server_response.read()
output_file = output_file_bytes.decode(encoding="UTF-8")


# In[ ]:


energy = []
len_energy = len("total energy")
for line in output_file.split("\n"):
    if line.strip().lstrip("!")[:len_energy] == "total energy":
        energy.append(float(line.split("=")[1].rstrip("Ry")))


# In[ ]:


# plot energy with iteration step
import matplotlib.pyplot as plt

get_ipython().run_line_magic('matplotlib', 'inline')

plt.plot(energy)
plt.xlabel("Number of iteration")
plt.ylabel("Energy (Ry)")
plt.show()


# In[ ]:




