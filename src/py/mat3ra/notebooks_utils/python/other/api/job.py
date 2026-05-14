import datetime
from collections import Counter
from typing import List

from mat3ra.api_client import JobEndpoints
from mat3ra.utils.extra.tabulate import pretty_print

from src.py.mat3ra.notebooks_utils.jupyterlite.interrupts import interruptible_polling_loop


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
        counts.get("submitted", 0) + counts.get("queued", 0),
        counts.get("active", 0),
        counts.get("finished", 0),
        counts.get("error", 0),
    ]
    pretty_print([row], headers, tablefmt="grid", stralign="center")

    active_statuses = {"pre-submission", "submitted", "queued", "active"}
    return not statuses or any(status in active_statuses for status in statuses)
