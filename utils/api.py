import datetime
import json
import os
import time
import urllib.request
from typing import List

from exabyte_api_client.endpoints.bank_workflows import BankWorkflowEndpoints
from exabyte_api_client.endpoints.jobs import JobEndpoints
from exabyte_api_client.endpoints.properties import PropertiesEndpoints
from tabulate import tabulate


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


def wait_for_jobs_to_finish(endpoint: JobEndpoints, job_ids: list, poll_interval: int = 10) -> None:
    """
    Waits for jobs to finish and prints their statuses.
    A job is considered finished if it is not in "pre-submission", "submitted", or, "active" status.

    Args:
        endpoint (JobEndpoints): Job endpoint object from the Exabyte API Client
        job_ids (list): list of job IDs to wait for
        poll_interval (int): poll interval for job information in seconds. Defaults to 10.
    """
    print("Wait for jobs to finish, poll interval: {0} sec".format(poll_interval))
    while True:
        statuses = get_jobs_statuses_by_ids(endpoint, job_ids)

        errored_jobs = len([status for status in statuses if status == "error"])
        active_jobs = len([status for status in statuses if status == "active"])
        finished_jobs = len([status for status in statuses if status == "finished"])
        submitted_jobs = len([status for status in statuses if status == "submitted"])

        headers = ["TIME", "SUBMITTED-JOBS", "ACTIVE-JOBS", "FINISHED-JOBS", "ERRORED-JOBS"]
        now = datetime.datetime.now().strftime("%Y-%m-%d-%H:%M:%S")
        row = [now, submitted_jobs, active_jobs, finished_jobs, errored_jobs]
        print(tabulate([row], headers, tablefmt="grid", stralign="center"))

        if all([status not in ["pre-submission", "submitted", "active"] for status in statuses]):
            break
        time.sleep(poll_interval)


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
