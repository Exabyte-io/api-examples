import datetime
import json
import os
import urllib.request
from collections import Counter
from typing import List, Optional, Union

from mat3ra.api_client import APIClient
from mat3ra.api_client.endpoints.bank_workflows import BankWorkflowEndpoints
from mat3ra.api_client.endpoints.jobs import JobEndpoints
from mat3ra.api_client.endpoints.properties import PropertiesEndpoints
from mat3ra.made.material import Material
from mat3ra.prode import PropertyName
from mat3ra.utils.extra.tabulate import pretty_print
from mat3ra.utils.jupyterlite.interrupts import interruptible_polling_loop
from mat3ra.wode import Workflow


def save_files(job_id: str, job_endpoint: JobEndpoints, filename_on_cloud: str, filename_on_disk: str) -> None:
    """
    Saves a file to disk, overwriting any files with the same name as filename_on_disk

    Args:
        job_id (str): ID of the job
        job_endpoint (JobEndpoints): Job endpoint object from the Exabyte API Client
        filename_on_cloud (str): Name of the file on the server
        filename_on_disk (str): Name the file will be saved to

    Returns:
        None
    """
    files = job_endpoint.list_files(job_id)
    for file in files:
        if file["name"] == filename_on_cloud:
            file_metadata = file

    # Get a download URL for the CONTCAR
    signed_url = file_metadata["signedUrl"]

    # Download the contcar to memory
    server_response = urllib.request.urlopen(signed_url)

    # Write it to disk
    with open(filename_on_disk, "wb") as outp:
        outp.write(server_response.read())


def get_jobs_statuses_by_ids(endpoint: JobEndpoints, job_ids: List[str]) -> List[str]:
    """
    Gets jobs statues by their IDs.

    Args:
        endpoint (JobEndpoints): Job endpoint object from the Exabyte API Client
        job_ids (list): list of job IDs to get the status for

    Returns:
        list: list of job statuses
    """
    jobs = endpoint.list({"_id": {"$in": job_ids}}, {"fields": {"status": 1}})
    return [job["status"] for job in jobs]


@interruptible_polling_loop()
def wait_for_jobs_to_finish_async(endpoint: JobEndpoints, job_ids: List[str]) -> bool:
    """
    Waits for jobs to finish and prints their statuses.
    A job is considered finished if it is not in "pre-submission", "submitted", or, "active" status.

    Args:
        endpoint (JobEndpoints): Job endpoint object from the Exabyte API Client
        job_ids (list): list of job IDs to wait for
    """
    statuses = get_jobs_statuses_by_ids(endpoint, job_ids)
    counts = Counter(statuses)
    headers = ["TIME", "SUBMITTED-JOBS", "ACTIVE-JOBS", "FINISHED-JOBS", "ERRORED-JOBS"]
    now = datetime.datetime.now().strftime("%Y-%m-%d-%H:%M:%S")
    row = [
        now,
        counts.get("submitted", 0),
        counts.get("active", 0),
        counts.get("finished", 0),
        counts.get("error", 0),
    ]
    pretty_print([row], headers, tablefmt="grid", stralign="center")

    active_statuses = {"pre-submission", "submitted", "active"}
    return any(status in active_statuses for status in statuses)


def copy_bank_workflow_by_system_name(endpoint: BankWorkflowEndpoints, system_name: str, account_id: str) -> dict:
    """
    Copies a bank workflow with given ID into the account's workflows.

    Args:
        endpoint (endpoints.bank_workflows.BankWorkflowEndpoints): an instance of BankWorkflowEndpoints class
        system_name (str): workflow system name.
        account_id (str): ID of account to copy the bank workflow into.

    Returns:
        dict: new account's workflow
    """
    bank_workflow_id = endpoint.list({"systemName": system_name})[0]["_id"]
    return endpoint.copy(bank_workflow_id, account_id)["_id"]


def get_property_by_subworkflow_and_unit_indicies(
    endpoint: PropertiesEndpoints, property_name: str, job: dict, subworkflow_index: int, unit_index: int
) -> dict:
    """
    Returns the property extracted in the given unit of the job's subworkflow.

    Args:
        endpoint (endpoints.properties.PropertiesEndpoints): an instance of PropertiesEndpoints class.
        property_name (str): name of property to extract.
        job (dict): job config to extract the property from.
        subworkflow_index (int): index of subworkflow to extract the property from.
        unit_index (int): index of unit to extract the property from.

    Returns:
        dict: extracted property
    """
    unit_flowchart_id = job["workflow"]["subworkflows"][subworkflow_index]["units"][unit_index]["flowchartId"]
    return endpoint.get_property(job["_id"], unit_flowchart_id, property_name)


def get_cluster_name(name: str = "cluster-001") -> str:
    clusters = json.loads(os.environ.get("CLUSTERS", "[]") or "[]")
    return clusters[0] if clusters else name


def get_or_create_material(api_client: APIClient, material, owner_id: str) -> dict:
    """
    Returns an existing material from the collection if one with the same structural hash
    exists under the given owner, otherwise creates a new one.
    Uses the client-side hash (mat3ra-made Material.hash) to avoid unnecessary DB writes.

    Args:
        api_client (APIClient): API client instance carrying the authorization context.
        material: mat3ra-made Material object (must have a .hash property).
        owner_id (str): Account ID under which to search and create.

    Returns:
        dict: The material dict (existing or newly created).
    """
    existing = api_client.materials.list({"hash": material.hash, "owner._id": owner_id})
    if existing:
        print(f"♻️  Reusing already existing Material: {existing[0]['_id']}")
        return existing[0]
    created = api_client.materials.create(material.to_dict(), owner_id=owner_id)
    print(f"✅ Material created: {created['_id']}")
    return created


def get_or_create_workflow(api_client: APIClient, workflow: Workflow, owner_id: str) -> dict:
    """
    Creates a workflow in the collection if none with the same hash exists under the given owner.
    Unit-level context (important settings) is stripped before saving so the base workflow
    stays clean and reusable.

    Args:
        api_client (APIClient): API client instance carrying the authorization context.
        workflow: mat3ra-wode Workflow object.
        owner_id (str): Account ID under which to search and create.

    Returns:
        dict: The workflow dict (existing or newly created).
    """
    existing = api_client.workflows.list({"hash": workflow.hash, "owner._id": owner_id})
    if existing:
        print(f"♻️  Reusing already existing Workflow: {existing[0]['_id']}")
        return existing[0]
    created = api_client.workflows.create(workflow.to_clean_dict(), owner_id=owner_id)
    print(f"✅ Workflow created: {created['_id']}")
    return created


FERMI_ENERGY_PROPERTIES = {
    PropertyName.non_scalar.band_structure.value,
    PropertyName.non_scalar.density_of_states.value,
}


def get_fermi_energy_for_job(client: APIClient, job_id: str, job: dict) -> Optional[float]:
    """
    Find the fermi_energy value for a job, mirroring calculateFermiEnergy in the web app.
    Finds the last unit in any subworkflow that declares fermi_energy in its results and returns its value.
    """
    for swf in job.get("workflow", {}).get("subworkflows", []):
        units_with_fe = [
            u
            for u in swf.get("units", [])
            if any(r.get("name") == PropertyName.scalar.fermi_energy.value for r in u.get("results", []))
        ]
        if not units_with_fe:
            continue
        last_unit_flowchart_id = units_with_fe[-1].get("flowchartId")
        fe_props = client.properties.get_for_job(job_id, PropertyName.scalar.fermi_energy.value, last_unit_flowchart_id)
        if fe_props:
            return fe_props[0].get("value")
    return None


def get_properties_for_job(client: APIClient, job_id: str, property_name: Optional[str] = None) -> List[dict]:
    """
    Fetch properties for a job, automatically enriching band_structure/DOS results with fermiEnergy.
    Use instead of client.properties.get_for_job when passing results to visualize_properties.
    """
    job = client.jobs.get(job_id)
    properties = client.properties.get_for_job(job_id, property_name)
    if property_name not in FERMI_ENERGY_PROPERTIES:
        return properties
    fermi_energy = get_fermi_energy_for_job(client, job_id, job)
    return [{**prop, "fermiEnergy": fermi_energy} for prop in properties]


def create_job(
    api_client: APIClient,
    materials: List[Union[dict, Material]],
    workflow: Union[dict, Workflow],
    project_id: str,
    owner_id: str,
    prefix: str,
    compute: Optional[dict] = None,
) -> List[dict]:
    """
    Creates jobs for each material using an embedded workflow with any context (important settings)
    already applied. The workflow _id is stripped so the server uses the embedded dict as-is,
    preserving unit-level context (kpath, kgrid, cutoffs, etc.) without saving them to the
    workflow collection.

    Args:
        api_client (APIClient): API client instance carrying the authorization context.
        materials (list): List of material dicts or mat3ra-made Material objects.
        workflow: Workflow dict or Workflow object with important settings already applied.
        project_id (str): Project ID.
        owner_id (str): Account ID.
        prefix (str): Job name prefix.
        compute (dict, optional): Compute configuration dict.

    Returns:
        list[dict]: List of created job dicts.
    """
    material_dicts = []
    for material in materials:
        if isinstance(material, Material):
            material_dicts.append(material.to_dict())
        else:
            material_dicts.append(material)

    job_workflow_dict = workflow.to_dict() if isinstance(workflow, Workflow) else workflow
    # Strip _id so the server uses the embedded workflow as-is instead of fetching from DB,
    # which would discard any unit-level context (kpath, kgrid, cutoffs, etc.).
    job_workflow_dict.pop("_id", None)
    is_multimaterial = job_workflow_dict.get("isMultimaterial", False)

    config = {
        "_project": {"_id": project_id},
        "workflow": job_workflow_dict,
        "owner": {"_id": owner_id},
        "name": prefix,
    }

    if is_multimaterial:
        config["_materials"] = [{"_id": mid} for mid in {md["_id"] for md in material_dicts}]
    else:
        config["_material"] = {"_id": material_dicts[0]["_id"]}

    if compute:
        config["compute"] = compute
    return api_client.jobs.create(config)
