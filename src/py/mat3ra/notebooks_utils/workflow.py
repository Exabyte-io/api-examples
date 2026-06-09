import re
from typing import Dict, List, Optional

from mat3ra.wode import Workflow


def patch_workflow_qe_input(
    workflow: Workflow,
    parameters: Dict[str, Dict[str, object]],
    unit_names: List[str],
    input_name: Optional[str] = None,
):
    """
    Patch QE inputs across workflow subworkflows for named units.

    Args:
        workflow: Workflow with subworkflows.
        parameters: Multi-section parameters {section: {key: val}}.
        unit_names: List of unit names to patch.
        input_name: Optional input file name filter.

    Example:
        patch_workflow_qe_input(workflow, {"system": {"vdw_corr": "d3_grimme"}}, ["pw_relax"])
    """
    f90 = lambda value: (  # noqa: E731
        f".{str(value).lower()}." if isinstance(value, bool) else repr(value) if isinstance(value, str) else str(value)
    )
    for subworkflow in workflow.subworkflows:
        for unit_name in unit_names:
            if not (unit := subworkflow.get_unit_by_name(name=unit_name)):
                continue
            for input_item in getattr(unit, "input", []):
                template = input_item.template
                if input_name not in (None, template.name):
                    continue
                content = template.content
                for section, updates in parameters.items():
                    name = section.lstrip("&")
                    match = re.search(rf"(?ims)(^&{re.escape(name)}\s*\n)(.*?)(^/\s*$)", content)
                    if not match:
                        raise ValueError(f"Namelist '&{name.upper()}' not found.")
                    header, body, footer = match.groups()
                    for key, value in updates.items():
                        line, pattern = f"    {key} = {f90(value)}", rf"(?im)^\s*{re.escape(key)}\s*=.*$"
                        body = re.sub(pattern, line, body) if re.search(pattern, body) else f"{body.rstrip()}\n{line}\n"
                    content = content[: match.start()] + header + body + footer + content[match.end() :]
                template.set_content(content)
            subworkflow.set_unit(unit)
