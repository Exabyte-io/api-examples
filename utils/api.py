import os
from typing import Any, Dict, List, Optional

import requests
from exabyte_api_client.endpoints.jobs import JobEndpoints
from exabyte_api_client.endpoints.materials import MaterialEndpoints
from exabyte_api_client.endpoints.projects import ProjectEndpoints
from exabyte_api_client.endpoints.workflows import WorkflowEndpoints

from .auth import ACCESS_TOKEN_ENV_VAR
from .settings import ACCOUNT_ID, ENDPOINT_ARGS

# TODO: remove in production
API_HOST = "localhost"
API_PORT = "3000"
API_SECURE = "false"
API_VERSION = "2018-10-01"


class Endpoints:
    def __init__(self, material, workflow, job, project):
        self.material = material
        self.workflow = workflow
        self.job = job
        self.project = project


def _create_endpoints_with_auth():
    # Create endpoints with ENDPOINT_ARGS
    material_endpoints = MaterialEndpoints(*ENDPOINT_ARGS)
    workflow_endpoints = WorkflowEndpoints(*ENDPOINT_ARGS)
    job_endpoints = JobEndpoints(*ENDPOINT_ARGS)
    project_endpoints = ProjectEndpoints(*ENDPOINT_ARGS)

    # Check if OIDC token is available
    access_token = os.environ.get(ACCESS_TOKEN_ENV_VAR)
    if access_token:
        # Use OIDC Bearer token authentication
        auth_header = {"Authorization": f"Bearer {access_token}"}
        material_endpoints.headers.update(auth_header)
        workflow_endpoints.headers.update(auth_header)
        job_endpoints.headers.update(auth_header)
        project_endpoints.headers.update(auth_header)

    return Endpoints(material_endpoints, workflow_endpoints, job_endpoints, project_endpoints)


endpoints = _create_endpoints_with_auth()


def get_owner_id() -> str:
    # If OIDC token exists, fetch real account ID
    access_token = os.environ.get(ACCESS_TOKEN_ENV_VAR)
    if access_token:
        protocol = "https" if API_SECURE else "http"
        port_str = f":{API_PORT}" if API_PORT not in [80, 443] else ""
        url = f"{protocol}://{API_HOST}{port_str}/api/v1/users/me"
        # Make request to /api/v1/users/me
        response = requests.get(url, headers={"Authorization": f"Bearer {access_token}"})
        account_id = response.json()["data"]["user"]["entity"]["defaultAccountId"]
        os.environ["ACCOUNT_ID"] = account_id  # Cache it
        return account_id  # Returns "k4WyvN62fTJPhzSpT" (real ID!)

    # Fallback to settings and validate
    account_id = os.environ.get("ACCOUNT_ID") or ACCOUNT_ID
    if account_id == "ACCOUNT_ID":
        raise ValueError("ACCOUNT_ID is not set. Please authenticate first")
    return account_id


def get_material(material_id: str) -> Dict[str, Any]:
    return endpoints.material.get(material_id)


def list_materials(query: Optional[Dict[str, Any]] = None, owner_id: Optional[str] = None) -> List[Dict[str, Any]]:
    owner = owner_id or get_owner_id()
    query_params = query or {}
    if "owner._id" not in query_params:
        query_params["owner._id"] = owner

    return endpoints.material.list(query_params)


def create_material(material: Any, owner_id: Optional[str] = None) -> Dict[str, Any]:
    owner = owner_id or get_owner_id()
    raw_config = material.to_dict()
    fields = ["name", "lattice", "basis"]
    config = {key: raw_config[key] for key in fields}

    return endpoints.material.create(config, owner)


def get_workflow(workflow_id: str) -> Dict[str, Any]:
    return endpoints.workflow.get(workflow_id)


def list_workflows(query: Optional[Dict[str, Any]] = None, owner_id: Optional[str] = None) -> List[Dict[str, Any]]:
    owner = owner_id or get_owner_id()
    query_params = query or {}
    if "owner._id" not in query_params:
        query_params["owner._id"] = owner

    return endpoints.workflow.list(query_params)


def create_workflow(workflow: Any, owner_id: Optional[str] = None) -> Dict[str, Any]:
    owner = owner_id or get_owner_id()
    config = workflow.to_dict()

    return endpoints.workflow.create(config, owner)


def get_job(job_id: str) -> Dict[str, Any]:
    return endpoints.job.get(job_id)


def list_jobs(query: Optional[Dict[str, Any]] = None, owner_id: Optional[str] = None) -> List[Dict[str, Any]]:
    owner = owner_id or get_owner_id()
    query_params = query or {}
    if "owner._id" not in query_params:
        query_params["owner._id"] = owner

    return endpoints.job.list(query_params)


def create_job(
    materials,
    workflow_id,
    project_id,
    name,
    compute,
    owner_id=None,
):
    owner = owner_id or get_owner_id()

    return endpoints.job.create_by_ids(
        materials=materials,
        workflow_id=workflow_id,
        project_id=project_id,
        owner_id=owner,
        prefix=name,
        compute=compute,
    )


def get_compute_config(cluster: str = "cluster-001") -> Dict[str, Any]:
    return endpoints.job.get_compute(cluster=cluster)


def submit_job(job_id: str) -> Dict[str, Any]:
    return endpoints.job.submit(job_id)


def get_default_project(owner_id: Optional[str] = None) -> str:
    owner = owner_id or get_owner_id()
    projects = endpoints.project.list({"isDefault": True, "owner._id": owner})

    if not projects:
        raise Exception(f"No default project found for owner {owner}")

    return projects[0]["_id"]
