from functools import wraps
from typing import Any, Callable, List

from mat3ra.api_client import APIClient
from mat3ra.made.material import Material
from mat3ra.wode import Workflow

from .core.entity.job.api import create_job as _create_job


# TODO: place to mat3ra-made
def convert_material_args(fn: Callable) -> Callable:
    """Converts any Material or List[Material] arg/kwarg to dict(s)."""

    @wraps(fn)
    def wrapper(*args: Any, **kwargs: Any) -> Any:
        def convert(v: Any) -> Any:
            if Material is not None and isinstance(v, Material):
                return v.to_dict()
            if isinstance(v, list):
                return [i.to_dict() if (Material is not None and isinstance(i, Material)) else i for i in v]
            return v

        return fn(*[convert(a) for a in args], **{k: convert(v) for k, v in kwargs.items()})

    return wrapper


# TODO: place to mat3ra-wode
def convert_workflow_args(fn: Callable) -> Callable:
    """Converts any Workflow arg/kwarg to a dict."""

    @wraps(fn)
    def wrapper(*args: Any, **kwargs: Any) -> Any:
        convert = lambda v: v.to_dict() if (Workflow is not None and isinstance(v, Workflow)) else v  # noqa: E731
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
