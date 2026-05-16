import datetime
from collections import Counter
from functools import wraps
from typing import Any, Callable, List

from mat3ra.api_client import APIClient, JobEndpoints
from mat3ra.made.material import Material
from mat3ra.utils.extra.tabulate import pretty_print
from mat3ra.wode import Workflow

from .core.entity.job.api import create_job as _create_job
from .core.entity.job.api import get_jobs_statuses_by_ids, save_files, submit_jobs
from .pyodide.runtime import interruptible_polling_loop


# TODO: place to mat3ra-made
def convert_material_args(fn: Callable) -> Callable:
    """Converts any Material or List[Material] arg/kwarg to dict(s)."""

    @wraps(fn)
    def wrapper(*args: Any, **kwargs: Any) -> Any:
        def convert(v: Any) -> Any:
            if isinstance(v, Material):
                return v.to_dict()
            if isinstance(v, list):
                return [i.to_dict() if isinstance(i, Material) else i for i in v]
            return v

        return fn(*[convert(a) for a in args], **{k: convert(v) for k, v in kwargs.items()})

    return wrapper


# TODO: place to mat3ra-wode
def convert_workflow_args(fn: Callable) -> Callable:
    """Converts any Workflow arg/kwarg to a dict."""

    @wraps(fn)
    def wrapper(*args: Any, **kwargs: Any) -> Any:
        convert = lambda v: v.to_dict() if isinstance(v, Workflow) else v  # noqa: E731
        return fn(*[convert(a) for a in args], **{k: convert(v) for k, v in kwargs.items()})

    return wrapper


create_job = convert_workflow_args(convert_material_args(_create_job))


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


__all__ = [
    "create_job",
    "get_convergence_series",
    "get_jobs_statuses_by_ids",
    "save_files",
    "submit_jobs",
    "wait_for_jobs_to_finish_async",
]
