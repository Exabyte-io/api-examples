from typing import List, Optional

from mat3ra.api_client import APIClient, PropertiesEndpoints
from mat3ra.prode import PropertyName

from .job import get_fermi_energy_flowchart_id

FERMI_ENERGY_PROPERTIES = {
    PropertyName.non_scalar.band_structure.value,
    PropertyName.non_scalar.density_of_states.value,
}


def get_properties_for_job(client: APIClient, job_id: str, property_name: Optional[str] = None) -> List[dict]:
    """
    Fetch properties for a job, automatically enriching band_structure/DOS results with fermiEnergy.
    Use instead of client.properties.get_for_job when passing results to visualize_properties.
    """
    job = client.jobs.get(job_id)
    properties = client.properties.get_for_job(job_id, property_name)
    if property_name not in FERMI_ENERGY_PROPERTIES:
        return properties
    flowchart_id = get_fermi_energy_flowchart_id(job)
    fermi_energy = None
    if flowchart_id:
        fe_props = client.properties.get_for_job(job_id, PropertyName.scalar.fermi_energy.value, flowchart_id)
        if fe_props:
            fermi_energy = fe_props[0].get("value")
    return [{**prop, "fermiEnergy": fermi_energy} for prop in properties]


def get_property_holder_for_job(
    client: APIClient, job_id: str, property_name: str, unit_id: Optional[str] = None
) -> dict:
    """
    Fetch the first full property holder for a job/property pair.

    Args:
        client (APIClient): API client instance.
        job_id (str): Job ID.
        property_name (str): Property name.
        unit_id (str, optional): Unit flowchart ID.

    Returns:
        dict: Full property holder document.
    """
    query = {
        "source.info.jobId": job_id,
        "data.name": property_name,
    }
    if unit_id:
        query["source.info.unitId"] = unit_id
    holders = client.properties.list(query=query)
    if not holders:
        raise ValueError(f"Property '{property_name}' not found for job '{job_id}'")
    return holders[0]


def update_property_holder_value(client: APIClient, property_holder_id: str, value: float) -> dict:
    """
    Update a scalar property's data.value.

    Args:
        client (APIClient): API client instance.
        property_holder_id (str): Property holder ID.
        value (float): New scalar value.

    Returns:
        dict: Server response payload.
    """
    return client.properties.update(property_holder_id, {"$set": {"data.value": value}})


def get_property_by_subworkflow_and_unit_indicies(
    endpoint: PropertiesEndpoints, property_name: str, job: dict, subworkflow_index: int, unit_index: int
) -> dict:
    """
    Returns the property extracted in the given unit of the job's subworkflow.

    Args:
        endpoint (PropertiesEndpoints): an instance of PropertiesEndpoints class.
        property_name (str): name of property to extract.
        job (dict): job config to extract the property from.
        subworkflow_index (int): index of subworkflow to extract the property from.
        unit_index (int): index of unit to extract the property from.

    Returns:
        dict: extracted property
    """
    unit_flowchart_id = job["workflow"]["subworkflows"][subworkflow_index]["units"][unit_index]["flowchartId"]
    return endpoint.get_property(job["_id"], unit_flowchart_id, property_name)
