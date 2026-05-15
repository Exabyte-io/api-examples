import datetime
import urllib.request
from collections import Counter
from typing import List, Optional, Union

from mat3ra.api_client import APIClient, JobEndpoints
from mat3ra.made.material import Material
from mat3ra.utils.extra.tabulate import pretty_print
from mat3ra.wode import Workflow

from ...pyodide.runtime import interruptible_polling_loop


def save_files(job_id: str, job_endpoint: JobEndpoints, filename_on_cloud: str, filename_on_disk: str) -> None:
    """
    Saves a file to disk, overwriting any files with the same name as filename_on_disk.

    Args:
        job_id (str): ID of the job
        job_endpoint (JobEndpoints): Job endpoint object from the Exabyte API Client
        filename_on_cloud (str): Name of the file on the server
        filename_on_disk (str): Name the file will be saved to
    """
    files = job_endpoint.list_files(job_id)
    file_metadata = next(f for f in files if filename_on_cloud in f["key"])
    signed_url = file_metadata["signedUrl"]
    server_response = urllib.request.urlopen(signed_url)
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
    A job is considered finished if it is not in "pre-submission", "submitted", or "active" status.

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
        counts.get("submitted", 0) + counts.get("queued", 0),
        counts.get("active", 0),
        counts.get("finished", 0),
        counts.get("error", 0),
    ]
    pretty_print([row], headers, tablefmt="grid", stralign="center")

    active_statuses = {"pre-submission", "submitted", "queued", "active"}
    return not statuses or any(status in active_statuses for status in statuses)


def create_job(
    api_client: APIClient,
    materials: List[Union[dict, Material]],
    workflow: Union[dict, Workflow],
    project_id: str,
    owner_id: str,
    prefix: str,
    compute: Optional[dict] = None,
) -> Union[dict, List[dict]]:
    """
    Creates jobs for each material using an embedded workflow with any context already applied.

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
    material_dicts = [m.to_dict() if isinstance(m, Material) else m for m in materials]
    job_workflow_dict = workflow.to_dict() if isinstance(workflow, Workflow) else workflow
    job_workflow_dict.pop("_id", None)
    is_multimaterial = job_workflow_dict.get("isMultiMaterial", False)

    config = {
        "_project": {"_id": project_id},
        "workflow": job_workflow_dict,
        "owner": {"_id": owner_id},
        "name": prefix,
    }

    if is_multimaterial:
        config["_material"] = {"_id": material_dicts[0]["_id"]}
        config["_materials"] = [{"_id": m["_id"]} for m in material_dicts]
    else:
        config["_material"] = {"_id": material_dicts[0]["_id"]}

    if compute:
        config["compute"] = compute
    return api_client.jobs.create(config)


def get_convergence_series(client: APIClient, job_id: str, subworkflow_index: int = 0) -> List[dict]:
    """
    Returns the convergence series from a finished convergence job.

    Args:
        client: API client instance.
        job_id: ID of the finished convergence job.
        subworkflow_index: Index of the convergence subworkflow (default 0).

    Returns:
        List of dicts with keys "x", "parameter", "y".
    """
    finished_job = client.jobs.get(job_id)
    job_workflow = Workflow.create(finished_job["workflow"])
    subworkflow = job_workflow.subworkflows[subworkflow_index]
    return subworkflow.convergence_series(finished_job.get("scopeTrack"))


def submit_jobs(endpoint: JobEndpoints, job_ids: List[str]) -> None:
    """
    Submits jobs by IDs.

    Args:
        endpoint (JobEndpoints): Job endpoint object from the Exabyte API Client.
        job_ids (list[str]): Job IDs to submit.
    """
    for job_id in job_ids:
        endpoint.submit(job_id)
