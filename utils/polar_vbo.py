from copy import deepcopy

from mat3ra.wode.workflows import Workflow

POLAR_POSTPROCESS_SUBWORKFLOW_ID = "df9ddd647d8540ef6fd947fa"
POLAR_POSTPROCESS_FLOWCHART_ID = "87b5f54edde75c917c0d9389"
CALCULATE_VBO_SUBWORKFLOW_NAME = "Calculate VBO"
POLAR_POSTPROCESS_SUBWORKFLOW_NAME = "Python VBO Polar PostProcess"

POLAR_VBO_SCRIPT = """# ------------------------------------------------------------------ #
# Linear Fit of ESP for Polar Interface VBO Calculation              #
# ------------------------------------------------------------------ #
#
# Reference: Choudhary & Garrity, arXiv:2401.02021 (InterMat)        #
#
# For polar interfaces, ESP shows linear gradient in bulk regions
# due to internal electric field. We fit each slab region and use
# the average value of the fit as the ESP reference.
#
# VBO Calculation:
#   1. Fit interface profile over slab 1 region -> Va_interface
#   2. Fit interface profile over slab 2 region -> Vb_interface
#   3. Fit bulk left profile over slab 1 region -> Va_bulk_left
#   4. Fit bulk right profile over slab 2 region -> Vb_bulk_right
#   5. VBO = (delta V_interface) - (delta V_bulk)
#      where delta V_interface = Vb_interface - Va_interface
#            delta V_bulk = Vb_bulk_right - Va_bulk_left
#
# Input:
#   - profile_left, profile_right: ESP profiles for bulk materials
#   - profile_interface: ESP profile for interface structure
#
# Output: VBO (Valence Band Offset)
# ------------------------------------------------------------------ #
import json
from types import SimpleNamespace

import ase.io
import matplotlib
import numpy as np
from mat3ra.made.material import Material
from mat3ra.made.tools.convert import from_ase
from scipy.stats import linregress

matplotlib.use("Agg")
import matplotlib.pyplot as plt

pw_scf_output = "./pw_scf.out"
pw_scf_output_1 = "./pw_scf.out-1"
pw_scf_output_2 = "./pw_scf.out-2"

atoms = ase.io.read(pw_scf_output, format="espresso-out")
atoms_1 = ase.io.read(pw_scf_output_1, format="espresso-out")
atoms_2 = ase.io.read(pw_scf_output_2, format="espresso-out")

material = Material.create(from_ase(atoms))
material_1 = Material.create(from_ase(atoms_1))
material_2 = Material.create(from_ase(atoms_2))

material.to_cartesian()
material_1.to_cartesian()
material_2.to_cartesian()

coords = material.basis.coordinates.values
elements = material.basis.elements.values
z_elements = sorted(zip([c[2] for c in coords], elements))
n_left = len(material_1.basis.elements.values)

z_max_1 = z_elements[n_left - 1][0]
z_min_2 = z_elements[n_left][0]
z_min_1 = z_elements[0][0]
z_max_2 = z_elements[-1][0]

print(f"Detected Slab 1 (left) boundaries: z = {z_min_1:.3f} to {z_max_1:.3f} A")
print(f"Detected Slab 2 (right) boundaries: z = {z_min_2:.3f} to {z_max_2:.3f} A")

checkpoint_file = "./.mat3ra/checkpoint"
with open(checkpoint_file, "r") as f:
    checkpoint_data = json.load(f)
    profile_interface = SimpleNamespace(
        **checkpoint_data["scope"]["local"]["average-electrostatic-potential"]["average_potential_profile"]
    )
    profile_left = SimpleNamespace(
        **checkpoint_data["scope"]["local"]["average-electrostatic-potential-left"]["average_potential_profile"]
    )
    profile_right = SimpleNamespace(
        **checkpoint_data["scope"]["local"]["average-electrostatic-potential-right"]["average_potential_profile"]
    )

X = np.array(profile_interface.xDataArray)
Y = np.array(profile_interface.yDataSeries[1])
X_left = np.array(profile_left.xDataArray)
Y_left = np.array(profile_left.yDataSeries[1])
X_right = np.array(profile_right.xDataArray)
Y_right = np.array(profile_right.yDataSeries[1])


def get_region_indices(x_data, x_min, x_max):
    mask = (x_data >= x_min) & (x_data <= x_max)
    indices = np.where(mask)[0]
    if len(indices) == 0:
        return 0, len(x_data)
    return indices[0], indices[-1] + 1


def fit_and_average(x_data, y_data, start_idx, end_idx):
    x_region = x_data[start_idx:end_idx]
    y_region = y_data[start_idx:end_idx]

    if len(x_region) < 2:
        avg = float(np.mean(y_region)) if len(y_region) > 0 else 0.0
        return avg, 0.0, avg

    slope, intercept, _, _, _ = linregress(x_region, y_region)
    x_mid = (x_region[0] + x_region[-1]) / 2.0
    avg_value = slope * x_mid + intercept
    return float(avg_value), float(slope), float(intercept)


slab1_start, slab1_end = get_region_indices(X, z_min_1, z_max_1)
slab2_start, slab2_end = get_region_indices(X, z_min_2, z_max_2)

Va_interface, slope_a, intercept_a = fit_and_average(X, Y, slab1_start, slab1_end)
Vb_interface, slope_b, intercept_b = fit_and_average(X, Y, slab2_start, slab2_end)

slab1_start_left, slab1_end_left = get_region_indices(X_left, z_min_1, z_max_1)
slab2_start_right, slab2_end_right = get_region_indices(X_right, z_min_2, z_max_2)

Va_bulk_left, _, _ = fit_and_average(X_left, Y_left, slab1_start_left, slab1_end_left)
Vb_bulk_right, _, _ = fit_and_average(X_right, Y_right, slab2_start_right, slab2_end_right)

VBO = (Vb_interface - Va_interface) - (Vb_bulk_right - Va_bulk_left)

print(f"Interface ESP Slab 1 (Va_interface): {Va_interface:.3f} eV")
print(f"Interface ESP Slab 2 (Vb_interface): {Vb_interface:.3f} eV")
print(f"Bulk ESP Left (Va_bulk): {Va_bulk_left:.3f} eV")
print(f"Bulk ESP Right (Vb_bulk): {Vb_bulk_right:.3f} eV")
print(f"Interface delta V: {Vb_interface - Va_interface:.3f} eV")
print(f"Bulk delta V: {Vb_bulk_right - Va_bulk_left:.3f} eV")
print(f"Valence Band Offset (VBO): {VBO:.3f} eV")

plt.figure(figsize=(10, 6))
plt.plot(X, Y, label="Macroscopic Average Potential", linewidth=2)
plt.axvspan(z_min_1, z_max_1, color="red", alpha=0.2, label="Slab 1 Region")
plt.axvspan(z_min_2, z_max_2, color="blue", alpha=0.2, label="Slab 2 Region")

if slab1_end > slab1_start:
    x_fit1 = X[slab1_start:slab1_end]
    y_fit1 = slope_a * x_fit1 + intercept_a
    plt.plot(x_fit1, y_fit1, color="darkred", linestyle="--", linewidth=2, label="Fit Slab 1")

if slab2_end > slab2_start:
    x_fit2 = X[slab2_start:slab2_end]
    y_fit2 = slope_b * x_fit2 + intercept_b
    plt.plot(x_fit2, y_fit2, color="darkblue", linestyle="--", linewidth=2, label="Fit Slab 2")

plt.axhline(Va_interface, color="red", linestyle=":", linewidth=2, label=f"Avg ESP Slab 1 = {Va_interface:.3f} eV")
plt.axhline(Vb_interface, color="blue", linestyle=":", linewidth=2, label=f"Avg ESP Slab 2 = {Vb_interface:.3f} eV")

plt.xlabel("z-coordinate (A)", fontsize=12)
plt.ylabel("Macroscopic Average Potential (eV)", fontsize=12)
plt.title(f"Polar Interface VBO = {VBO:.3f} eV", fontsize=14, fontweight="bold")
plt.legend(loc="best", fontsize=10)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("polar_vbo_fit_interface.png", dpi=150, bbox_inches="tight")
plt.close()
"""

POLAR_VBO_REQUIREMENTS = """munch==2.5.0
numpy>=1.19.5
scipy>=1.5.4
matplotlib>=3.0.0
ase>=3.22.1
mat3ra-made[tools]>=2024.11.12.post0
"""


def _polar_postprocess_subworkflow(source_python_subworkflow):
    source_unit = next(unit for unit in source_python_subworkflow["units"] if unit["type"] == "execution")
    source_flavor_inputs = source_unit.get("flavor", {}).get("input", [])
    script_input_name = source_flavor_inputs[0].get("name", "script.py") if source_flavor_inputs else "script.py"
    script_template_name = (
        source_flavor_inputs[0].get("templateName", script_input_name) if source_flavor_inputs else script_input_name
    )
    requirements_input_name = (
        source_flavor_inputs[1].get("name", "requirements.txt") if len(source_flavor_inputs) > 1 else "requirements.txt"
    )
    requirements_template_name = (
        source_flavor_inputs[1].get("templateName", requirements_input_name)
        if len(source_flavor_inputs) > 1
        else requirements_input_name
    )

    polar_unit = deepcopy(source_unit)
    polar_unit["flowchartId"] = "20ea6fb714ad8529fae09fd3"
    polar_unit["name"] = "Get VBO for Polar material"
    polar_unit["head"] = True
    polar_unit.pop("next", None)
    polar_unit["context"] = {"subworkflowContext": {}}
    polar_unit["results"] = [{"name": "file_content", "filetype": "image", "basename": "polar_vbo_fit_interface.png"}]
    polar_unit["input"] = [
        {
            "content": POLAR_VBO_SCRIPT,
            "_id": "cKsDNiq3STXgGX9hP",
            "applicationName": "python",
            "executableName": "python",
            "name": script_input_name,
            "templateName": script_template_name,
            "contextProviders": [],
            "rendered": POLAR_VBO_SCRIPT,
            "schemaVersion": "2022.8.16",
        },
        {
            "content": POLAR_VBO_REQUIREMENTS,
            "_id": "yTM9FqBpxmFf2iigk",
            "applicationName": "python",
            "executableName": "python",
            "name": requirements_input_name,
            "templateName": requirements_template_name,
            "contextProviders": [],
            "rendered": POLAR_VBO_REQUIREMENTS,
            "schemaVersion": "2022.8.16",
        },
    ]

    polar_subworkflow = deepcopy(source_python_subworkflow)
    polar_subworkflow["_id"] = POLAR_POSTPROCESS_SUBWORKFLOW_ID
    polar_subworkflow["name"] = POLAR_POSTPROCESS_SUBWORKFLOW_NAME
    polar_subworkflow["properties"] = []
    polar_subworkflow["units"] = [polar_unit]
    return polar_subworkflow


def add_polar_vbo_postprocess(workflow: Workflow) -> Workflow:
    workflow_dict = deepcopy(workflow.to_dict())

    if any(unit["name"] == POLAR_POSTPROCESS_SUBWORKFLOW_NAME for unit in workflow_dict["units"]):
        return workflow

    calculate_vbo_unit = next(unit for unit in workflow_dict["units"] if unit["name"] == CALCULATE_VBO_SUBWORKFLOW_NAME)
    calculate_vbo_unit["next"] = POLAR_POSTPROCESS_FLOWCHART_ID

    source_python_subworkflow = next(
        subworkflow for subworkflow in workflow_dict["subworkflows"] if subworkflow["application"]["name"] == "python"
    )

    workflow_dict["units"].append(
        {
            "name": POLAR_POSTPROCESS_SUBWORKFLOW_NAME,
            "type": "subworkflow",
            "_id": POLAR_POSTPROCESS_SUBWORKFLOW_ID,
            "flowchartId": POLAR_POSTPROCESS_FLOWCHART_ID,
            "status": "idle",
            "statusTrack": [],
            "tags": [],
            "head": False,
        }
    )
    workflow_dict["subworkflows"].append(_polar_postprocess_subworkflow(source_python_subworkflow))
    return Workflow.create(workflow_dict)
