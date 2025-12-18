import os
from typing import Any, Dict, List, Optional

from exabyte_api_client.endpoints.jobs import JobEndpoints
from exabyte_api_client.endpoints.materials import MaterialEndpoints
from exabyte_api_client.endpoints.projects import ProjectEndpoints
from exabyte_api_client.endpoints.workflows import WorkflowEndpoints

from .auth import ACCESS_TOKEN_ENV_VAR
from .settings import ACCOUNT_ID, ENDPOINT_ARGS

API_HOST = "localhost"
API_PORT = "3000"
API_SECURE = "false"
API_VERSION = "2018-10-01"


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

    return material_endpoints, workflow_endpoints, job_endpoints, project_endpoints


def get_owner_id() -> str:
    import requests

    access_token = os.environ.get(ACCESS_TOKEN_ENV_VAR)

    if access_token:
        # Use OIDC to get account ID from /users/me endpoint
        try:
            from .settings import HOST, PORT, SECURE

            protocol = "https" if SECURE else "http"
            port_str = f":{PORT}" if PORT not in [80, 443] else ""
            url = f"{protocol}://{HOST}{port_str}/api/v1/users/me"

            headers = {"Authorization": f"Bearer {access_token}", "Content-Type": "application/json"}

            response = requests.get(url, headers=headers, timeout=30)
            response.raise_for_status()
            user_data = response.json()

            # Extract account ID from response
            account_id = user_data["data"]["user"]["entity"]["defaultAccountId"]

            # Cache it in environment for future use
            os.environ["ACCOUNT_ID"] = account_id

            return account_id
        except Exception as e:
            # If OIDC fetch fails, fall back to settings
            print(f"⚠️  Failed to get account ID from OIDC: {e}")
            print("⚠️  Falling back to ACCOUNT_ID from settings")

    # Return ACCOUNT_ID from settings (works for traditional auth)
    account_id = os.environ.get("ACCOUNT_ID") or ACCOUNT_ID

    # Validate it's not the placeholder string
    if account_id == "ACCOUNT_ID":
        raise ValueError(
            "ACCOUNT_ID is not set. Please authenticate first using: await authenticate() "
            "or set ACCOUNT_ID in your settings.json file."
        )

    return account_id


def get_material(material_id: str) -> Dict[str, Any]:
    material_endpoints, _, _, _ = _create_endpoints_with_auth()
    return material_endpoints.get(material_id)


def list_materials(query: Optional[Dict[str, Any]] = None, owner_id: Optional[str] = None) -> List[Dict[str, Any]]:
    material_endpoints, _, _, _ = _create_endpoints_with_auth()
    owner = owner_id or get_owner_id()

    query_params = query or {}
    if "owner._id" not in query_params:
        query_params["owner._id"] = owner

    return material_endpoints.list(query_params)


def create_material(material: Any, owner_id: Optional[str] = None) -> Dict[str, Any]:
    material_endpoints, _, _, _ = _create_endpoints_with_auth()
    owner = owner_id or get_owner_id()

    raw_config = material.to_dict()
    fields = ["name", "lattice", "basis"]
    config = {key: raw_config[key] for key in fields}

    return material_endpoints.create(config, owner)


def get_workflow(workflow_id: str) -> Dict[str, Any]:
    _, workflow_endpoints, _, _ = _create_endpoints_with_auth()
    return workflow_endpoints.get(workflow_id)


def list_workflows(query: Optional[Dict[str, Any]] = None, owner_id: Optional[str] = None) -> List[Dict[str, Any]]:
    _, workflow_endpoints, _, _ = _create_endpoints_with_auth()
    owner = owner_id or get_owner_id()

    query_params = query or {}
    if "owner._id" not in query_params:
        query_params["owner._id"] = owner

    return workflow_endpoints.list(query_params)


def create_workflow(workflow: Any, owner_id: Optional[str] = None) -> Dict[str, Any]:
    _, workflow_endpoints, _, _ = _create_endpoints_with_auth()
    owner = owner_id or get_owner_id()

    config = workflow.to_dict()

    # Clean up None values that cause API errors
    config.pop("metadata", None)
    config.pop("compute", None)

    for swf in config.get("subworkflows", []):
        swf.pop("compute", None)
        if "model" in swf and "method" in swf["model"]:
            swf["model"]["method"].pop("precision", None)

    for unit in config.get("units", []):
        unit.pop("status", None)
        unit.pop("tags", None)
        unit.pop("statusTrack", None)

    return workflow_endpoints.create(config, owner)


def get_job(job_id: str) -> Dict[str, Any]:
    _, _, job_endpoints, _ = _create_endpoints_with_auth()
    return job_endpoints.get(job_id)


def list_jobs(query: Optional[Dict[str, Any]] = None, owner_id: Optional[str] = None) -> List[Dict[str, Any]]:
    _, _, job_endpoints, _ = _create_endpoints_with_auth()
    owner = owner_id or get_owner_id()

    query_params = query or {}
    if "owner._id" not in query_params:
        query_params["owner._id"] = owner

    return job_endpoints.list(query_params)


def create_job(
    materials,
    workflow_id,
    project_id,
    name,
    compute,
    owner_id=None,
):
    _, _, job_endpoints, _ = _create_endpoints_with_auth()
    owner = owner_id or get_owner_id()

    return job_endpoints.create_by_ids(
        materials=materials,
        workflow_id=workflow_id,
        project_id=project_id,
        owner_id=owner,
        prefix=name,
        compute=compute,
    )


def get_compute_config(cluster: str = "cluster-001") -> Dict[str, Any]:
    _, _, job_endpoints, _ = _create_endpoints_with_auth()
    return job_endpoints.get_compute(cluster=cluster)


def submit_job(job_id: str) -> Dict[str, Any]:
    _, _, job_endpoints, _ = _create_endpoints_with_auth()
    return job_endpoints.submit(job_id)


def get_default_project(owner_id: Optional[str] = None) -> str:
    _, _, _, project_endpoints = _create_endpoints_with_auth()
    owner = owner_id or get_owner_id()

    projects = project_endpoints.list({"isDefault": True, "owner._id": owner})

    if not projects:
        raise Exception(f"No default project found for owner {owner}")

    return projects[0]["_id"]
