#!/usr/bin/env python
# coding: utf-8

# <a href="https://colab.research.google.com/github/Exabyte-io/api-examples/blob/dev/examples/workflow/qe_scf_calculation.ipynb" target="_parent">
# <img alt="Open in Google Colab" src="https://user-images.githubusercontent.com/20477508/128780728-491fea90-9b23-495f-a091-11681150db37.jpeg" width="150" border="0">
# </a>

# # Quantum Espresso SCF calcuation via API
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


# In[1]:


from utils.settings import ENDPOINT_ARGS, ACCOUNT_ID
from utils.generic import display_JSON

from exabyte_api_client.endpoints.workflows import WorkflowEndpoints


# In[2]:


# Initialize a helper class to interact with WorkflowEndpoints
endpoint = WorkflowEndpoints(*ENDPOINT_ARGS)


# In[3]:


# payload for workflow creation
BODY = {
    "name": "Silicon-SCF-bash-REST",
    "subworkflows": [
        {
            "name": "QE-SCF",
            "application": {
                "name": "shell",
                "version": "4.2.46",
                "build": "Default",
                "isDefault": "true",
                "summary": "Shell Script",
                "shortName": "sh",
            },
            "properties": ["total_energy", "fermi_energy"],
            "units": [
                {
                    "type": "execution",
                    "application": {
                        "name": "shell",
                        "version": "4.2.46",
                        "build": "Default",
                        "isDefault": "true",
                        "summary": "Shell Script",
                        "shortName": "sh",
                    },
                    "head": "true",
                    "input": [
                        {
                            "content": "#!/bin/bash\n# ---------------------- CLUSTER PARAMETERS ---------------------- #\n#PBS -N Silicon-SCF\n#PBS -j oe\n#PBS -l nodes=1\n#PBS -l ppn=4\n#PBS -q OR\n#PBS -l walltime=00:00:30:00\n#PBS -A seminar-espresso-tutorials\n# ------------------------- INPUT FILES -------------------------- #\n# switch to the job working directory\ncd $PBS_O_WORKDIR\n\nBASE_URL='https://github.com/pranabdas/espresso/raw/next/src'\nIN_FILE='pw.scf.silicon.in'\nPP_FILE='Si.pz-vbc.UPF'\n\n# get QE input file\nwget ${BASE_URL}/silicon/${IN_FILE} -O ${IN_FILE}\n\n# get pseudopotential file\nwget ${BASE_URL}/pseudos/${PP_FILE} -O ${PP_FILE}\n\n# delete pseudo_dir from input file and provide via env\nsed -i '/^\\s*pseudo_dir/d' ${IN_FILE}\nexport ESPRESSO_PSEUDO='./'\n\n# --------------------------- RUN JOB ---------------------------- #\n# load required module\nmodule add espresso/63-i-174-impi-044\n\n# run the calculation\nmpirun -np $PBS_NP pw.x -in ${IN_FILE} | tee ${IN_FILE%*.in}.out\n\n# ----------------------------- END ------------------------------ #\n",
                            "name": "pbs-Silicon-scf.sh",
                        }
                    ],
                    "name": "Silicon-SCF.sh",
                    "executable": {
                        "isDefault": "true",
                        "monitors": ["standard_output"],
                        "results": [
                            "atomic_forces",
                            "band_gaps",
                            "band_structure",
                            "density_of_states",
                            "fermi_energy",
                            "phonon_dispersions",
                            "phonon_dos",
                            "pressure",
                            "stress_tensor",
                            "total_energy",
                            "total_energy_contributions",
                            "total_force",
                            "zero_point_energy",
                            "final_structure",
                            "magnetic_moments",
                            "reaction_energy_barrier",
                            "reaction_energy_profile",
                            "potential_profile",
                            "charge_density_profile",
                        ],
                        "name": "sh",
                    },
                    "flavor": {
                        "isDefault": "true",
                        "input": [{"name": "hello_world.sh"}],
                        "monitors": ["standard_output"],
                        "applicationName": "shell",
                        "executableName": "sh",
                        "name": "hello_world",
                        "executable": {
                            "isDefault": "true",
                            "monitors": ["standard_output"],
                            "results": [
                                "atomic_forces",
                                "band_gaps",
                                "band_structure",
                                "density_of_states",
                                "fermi_energy",
                                "phonon_dispersions",
                                "phonon_dos",
                                "pressure",
                                "stress_tensor",
                                "total_energy",
                                "total_energy_contributions",
                                "total_force",
                                "zero_point_energy",
                                "final_structure",
                                "magnetic_moments",
                                "reaction_energy_barrier",
                                "reaction_energy_profile",
                                "potential_profile",
                                "charge_density_profile",
                            ],
                            "name": "sh",
                        },
                    },
                }
            ],
            "compute": {
                "ppn": 4,
                "nodes": 1,
                "queue": "OR",
                "timeLimit": "01:00:00",
                "notify": "n",
                "cluster": {"fqdn": "master-production-20160630-cluster-001.exabyte.io"},
            },
        }
    ],
}


# In[4]:


response = endpoint.create(BODY)


# In[ ]:




