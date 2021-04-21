#!/usr/bin/env python
# coding: utf-8

# # Overview
# 
# This example was created as part of our [Advanced Topics Webinar](https://www.youtube.com/watch?v=psSFC409jSg) on February 19, 2021. This webinar focused on explaining our API in detail, and provided examples of many areas of its functionality.
# 
# In this notebook, we showcase a major advantage of APIs: interoperability. We begin by performing a query using the [Materials Project](https://materialsproject.org) API for all systems containing Iron and Oxygen. We then filter our results (for demonstraiton purposes, we keep only the first 10 materials found). Finally, we upload our results to the Exabyte platform, where further calculations could be performed to characterize these materials.

# # Execution
# 
# > <span style="color: orange">**NOTE**</span>: In order to run this example, an active Exabyte.io account is required. RESTful API credentials shall be updated in [settings](../settings.py). The generation of the credentials is also explained therein.
# 
# ## Import packages

# In[]:


import os
import sys

# Import settings file and utils
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)
from utils.generic import ensure_packages_are_installed, display_JSON
ensure_packages_are_installed()
from settings import ENDPOINT_ARGS, MATERIALS_PROJECT_API_KEY

import ase.io
import pymatgen

from exabyte_api_client.endpoints.materials import MaterialEndpoints


# # Query the Materials Project
# 
# We begin by using the Materials Project API implemented in [PyMatGen](https://pymatgen.org/pymatgen.ext.matproj.html) to perform a query for all systems containing Iron and Oxygen.

# # Query Materials Project for all systems containing Iron and Oxygen
# materials_project_api = pymatgen.ext.matproj.MPRester(MATERIALS_PROJECT_API_KEY)
# iron_oxides_ids = materials_project_api.get_materials_ids("Fe-O")
# 
# print(iron_oxides_ids)

# # Filtering the Results
# 
# This returns a lot of materials - 160 to be exact! In many cases, it is useful to filter down the number of materials. For example, we may want to exclude large unit cells that may be computationally intensive to study. Or we may want to restrict our results to only thermodynamically-stable oxides, by use of the material's [energy above hull](https://wiki.materialsproject.org/Glossary_of_Terms).
# 
# As a basic example, here we only keep the first 10 iron oxides that the Materials Project API returned to us, and discard the other 150.

# In[]:


#As a demonstration, take the first 10 iron oxides
some_iron_oxides = iron_oxides_ids[:10]
print(some_iron_oxides)


# # Bringing Materials Into the User Account
# 
# Now that we have filtered the results from Materials Project down to just 10 structures, we may want to study them further with the computational models provided by Exabyte. For example, we may be interested in leveraging a DFT code to find the structure with the largest band-gap, or perhaps we want to conduct a high-throughput screening of each material's surface energies.

# In[]:


# Upload the first 10 iron oxides found to our account
exabyte_materials_api = MaterialEndpoints(*ENDPOINT_ARGS)
materials = exabyte_materials_api.import_from_materialsproject(MATERIALS_PROJECT_API_KEY, some_iron_oxides)


# Finally, it is always useful to stay organized. Materials sets make this convenient, acting as a folder to keep a group of related materials in. This would be especially helpful if, in the future, we wanted run a calculation over all the oxides we found in this example.

# In[]:


# Move the iron oxides to a materials set, just for this example
materials_set = exabyte_materials_api.create_set({"name" : "Some Iron Oxides"})
for material in materials:
    exabyte_materials_api.move_to_set(material["_id"], "", materials_set["_id"])

