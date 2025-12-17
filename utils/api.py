import os
from typing import Any, Dict, List, Optional, Union

import requests
from exabyte_api_client.endpoints.materials import MaterialEndpoints
from exabyte_api_client.endpoints.projects import ProjectEndpoints
from exabyte_api_client.endpoints.workflows import WorkflowEndpoints

from .auth import ACCESS_TOKEN_ENV_VAR
from .settings import ACCOUNT_ID, ENDPOINT_ARGS, HOST, PORT, SECURE, VERSION


def create_material(material: Any, owner_id: Optional[str] = None) -> Dict[str, Any]:
    """
    Creates a material on the platform.

    Args:
        material: Material object with a to_dict() method
        owner_id: Optional owner ID. Defaults to ACCOUNT_ID from settings

    Returns:
        dict: Created material response from the API
    """
    endpoint = MaterialEndpoints(*ENDPOINT_ARGS)
    owner = owner_id or ACCOUNT_ID

    raw_config = material.to_dict()
    fields = ["name", "lattice", "basis"]
    config = {key: raw_config[key] for key in fields}

    return endpoint.create(config, owner)


def create_workflow(workflow: Any, owner_id: Optional[str] = None) -> Dict[str, Any]:
    """
    Creates a workflow on the platform.

    Args:
        workflow: Workflow object with a to_dict() method
        owner_id: Optional owner ID. Defaults to ACCOUNT_ID from settings

    Returns:
        dict: Created workflow response from the API
    """
    endpoint = WorkflowEndpoints(*ENDPOINT_ARGS)
    owner = owner_id or ACCOUNT_ID

    config = workflow.to_dict()

    return endpoint.create(config, owner)


def get_default_project(owner_id: Optional[str] = None) -> str:
    """
    Gets the default project ID for an owner.

    Args:
        owner_id: Optional owner ID. Defaults to ACCOUNT_ID from settings

    Returns:
        str: Default project ID
    """
    endpoint = ProjectEndpoints(*ENDPOINT_ARGS)
    owner = owner_id or ACCOUNT_ID

    projects = endpoint.list({"isDefault": True, "owner._id": owner})

    return projects[0]["_id"]


def _get_api_base_url() -> str:
    protocol = "https" if SECURE else "http"
    port_str = f":{PORT}" if PORT not in [80, 443] else ""
    return f"{protocol}://{HOST}{port_str}/api/{VERSION}"


def _make_oidc_request(method: str, endpoint: str, data: Optional[Dict] = None) -> Union[Dict[str, Any], List[Any]]:
    access_token = os.environ.get(ACCESS_TOKEN_ENV_VAR)
    if not access_token:
        raise Exception(
            f"No OIDC access token found in environment variable '{ACCESS_TOKEN_ENV_VAR}'. "
            "Please authenticate first using: await authenticate()"
        )
    base_url = _get_api_base_url()
    url = f"{base_url}/{endpoint}"
    headers = {"Authorization": f"Bearer {access_token}", "Content-Type": "application/json"}
    response = requests.request(method=method, url=url, headers=headers, json=data, timeout=30)
    response.raise_for_status()
    return response.json()


def create_material_oidc(material: Any, owner_id: Optional[str] = None) -> Dict[str, Any]:
    owner = owner_id or ACCOUNT_ID
    raw_config = material.to_dict()
    fields = ["name", "lattice", "basis"]
    config = {key: raw_config[key] for key in fields}
    config["owner"] = {"_id": owner}
    return _make_oidc_request("PUT", "materials/create", config)


def create_workflow_oidc(workflow: Any, owner_id: Optional[str] = None) -> Dict[str, Any]:
    owner = owner_id or ACCOUNT_ID
    config = workflow.to_dict()
    config["owner"] = {"_id": owner}
    return _make_oidc_request("PUT", "workflows/create", config)


def get_default_project_oidc(owner_id: Optional[str] = None) -> str:
    owner = owner_id or ACCOUNT_ID
    response = _make_oidc_request("GET", f"projects?isDefault=true&owner._id={owner}")
    if isinstance(response, list):
        projects = response
    else:
        projects = response.get("data", [])
    return projects[0]["_id"]
