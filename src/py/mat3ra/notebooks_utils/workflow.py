import re
from typing import List, Mapping, Optional


def _format_to_f90_value(value: object) -> str:
    """Format Python value as Fortran namelist value."""
    if isinstance(value, bool):
        return ".true." if value else ".false."
    return f"'{value}'" if isinstance(value, str) else str(value)


def set_content(content: str, section: str, parameters: Mapping[str, object]) -> str:
    """Upsert parameters into a QE namelist section."""
    section_name = section.lstrip("&")
    pattern = rf"(?ims)(^&{re.escape(section_name)}\s*\n)(.*?)(^/\s*$)"
    match = re.search(pattern, content)
    if not match:
        raise ValueError(f"Namelist '&{section_name.upper()}' not found in input template.")

    before, header, body, footer, after = content[: match.start()], *match.groups(), content[match.end() :]

    for param, value in parameters.items():
        line = f"    {param} = {_format_to_f90_value(value)}"
        param_pattern = rf"(?im)^\s*{re.escape(param)}\s*=.*$"
        body = re.sub(param_pattern, line, body) if re.search(param_pattern, body) else body.rstrip() + f"\n{line}\n"

    return before + header + body + footer + after


def patch_qe_input(
    unit,
    parameters: Mapping[str, Mapping[str, object]],
    input_name: Optional[str] = None,
) -> None:
    """
    Patch QE namelist parameters on a workflow unit.

    Args:
        unit: Execution unit with input templates.
        parameters: Namelist parameters as {section: {key: value}}.
        input_name: Optional input file name filter.

    Example:
        patch_qe_input(unit, {"system": {"vdw_corr": "d3_grimme"}})
    """
    matched = False
    for item in getattr(unit, "input", []):
        template = item.template
        if input_name and template.name != input_name:
            continue

        content = template.content
        for section, params in parameters.items():
            content = set_content(content, section, params)
        template.set_content(content)
        matched = True

    if not matched:
        raise ValueError("No matching input template found for QE patch.")


def patch_workflow_qe_input(
    workflow,
    parameters: Mapping[str, Mapping[str, object]],
    unit_names: List[str],
    input_name: Optional[str] = None,
) -> None:
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
    for subworkflow in workflow.subworkflows:
        for unit_name in unit_names:
            if unit := subworkflow.get_unit_by_name(name=unit_name):
                patch_qe_input(unit, parameters, input_name=input_name)
                subworkflow.set_unit(unit)
