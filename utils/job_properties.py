from typing import Optional

from mat3ra.prode import PropertyName


def get_fermi_energy_flowchart_id(job: dict) -> Optional[str]:
    """
    Find the flowchart ID of the last unit that declares fermi_energy in its results.
    Mirrors calculateFermiEnergy logic in the web app.
    """
    for subworkflow in job.get("workflow", {}).get("subworkflows", []):
        units_with_fermi_energy = [
            unit
            for unit in subworkflow.get("units", [])
            if any(r.get("name") == PropertyName.scalar.fermi_energy.value for r in unit.get("results", []))
        ]
        if not units_with_fermi_energy:
            continue
        return units_with_fermi_energy[-1].get("flowchartId")
    return None
