from typing import Any, Dict, Optional, Union

from exabyte_api_client.endpoints.materials import MaterialEndpoints
from exabyte_api_client.endpoints.projects import ProjectEndpoints
from exabyte_api_client.endpoints.workflows import WorkflowEndpoints
from mat3ra.wode.workflows import Workflow

from .settings import ACCOUNT_ID, ENDPOINT_ARGS


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

    if isinstance(material, dict):
        raw_config = material
    else:
        raw_config = material.to_dict()
    fields = ["name", "lattice", "basis"]
    config = {key: raw_config[key] for key in fields}

    return endpoint.create(config, owner)


def create_workflow(workflow: Union[Workflow, Dict], owner_id: Optional[str] = None) -> Dict[str, Any]:
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
    if isinstance(workflow, dict):
        config = workflow
    else:
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
