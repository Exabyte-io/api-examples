import io
import re
from typing import Any, Dict, List, Mapping, Optional, Tuple

import f90nml
from mat3ra.utils.extra.jinja import JINJA_EXPRESSION_PATTERN


def _protect_jinja(text: str) -> Tuple[str, Dict[str, str]]:
    """Replace Jinja expressions with f90nml-safe placeholders."""
    placeholders: Dict[str, str] = {}
    pattern = re.compile(JINJA_EXPRESSION_PATTERN + r"|\{%-?.*?-?%\}", re.DOTALL)

    def replacer(match):
        key = f"__JINJA_{len(placeholders)}__"
        placeholders[key] = match.group(0)
        return f"'{key}'"

    return pattern.sub(replacer, text), placeholders


def _restore_jinja(text: str, placeholders: Mapping[str, str]) -> str:
    for key, original in placeholders.items():
        text = text.replace(f"'{key}'", original)
    return text


def set_content(content: str, section: str, parameters: Mapping[str, Any]) -> str:
    """Upsert parameters into a QE namelist section, preserving Jinja and case."""
    normalized = section.lstrip("&").upper()
    match = re.search(rf"(?ms)(^&{re.escape(normalized)}\s*\n.*?^/\s*$)", content)
    if not match:
        raise ValueError(f"Namelist '&{normalized}' not found in input template.")

    before, block, after = content[: match.start()], match.group(0), content[match.end() :]
    protected, placeholders = _protect_jinja(block)
    nml = f90nml.reads(protected)
    nml.patch({normalized.lower(): dict(parameters)})
    buffer = io.StringIO()
    nml.write(buffer)
    patched = re.sub(r"^&\w+", f"&{normalized}", buffer.getvalue(), count=1, flags=re.MULTILINE)
    return before + _restore_jinja(patched, placeholders) + after


def _get_template_attr(item, attr: str):
    """Get attribute from nested template dict/object or flat stub."""
    if isinstance(item, dict):
        template = item.get("template", item)
        return template.get(attr) if isinstance(template, dict) else None
    template = getattr(item, "template", item)
    return getattr(template, attr, None)


def _set_template_content(item, content: str):
    """Set content on nested template dict/object or flat stub."""
    if isinstance(item, dict):
        template = item.get("template", item)
        (template if isinstance(template, dict) else item)["content"] = content
    else:
        getattr(item, "template", item).content = content


def _patch_unit(unit, section: str, parameters: Mapping[str, Any], input_name: Optional[str]):
    """Patch a single namelist section across matching unit inputs."""
    matched = False
    for item in getattr(unit, "input", []):
        content = _get_template_attr(item, "content")
        if content and (not input_name or _get_template_attr(item, "name") == input_name):
            _set_template_content(item, set_content(content, section, parameters))
            matched = True
    if not matched:
        raise ValueError("No matching input template found for QE patch.")


def patch_qe_input(
    unit, section_or_params: Any, parameters: Optional[Mapping[str, Any]] = None, input_name: Optional[str] = None
):
    """
    Patch QE namelist parameters on a workflow unit.

    Examples:
        patch_qe_input(unit, "system", {"vdw_corr": "d3_grimme"})
        patch_qe_input(unit, {"system": {"vdw_corr": "d3_grimme"}})
    """
    if parameters is None:
        if not isinstance(section_or_params, Mapping) or not all(
            isinstance(v, Mapping) for v in section_or_params.values()
        ):
            raise TypeError("Expected section name with parameters or multi-section dict.")
        for section, params in section_or_params.items():
            _patch_unit(unit, section, params, input_name)
    else:
        _patch_unit(unit, section_or_params, parameters, input_name)


def patch_workflow_qe_input(
    workflow, parameters: Mapping[str, Mapping[str, Any]], unit_names: List[str], input_name: Optional[str] = None
):
    """
    Patch QE inputs across workflow subworkflows for named units.

    Example:
        patch_workflow_qe_input(workflow, {"system": {"vdw_corr": "d3_grimme"}}, ["pw_relax"])
    """
    for subworkflow in workflow.subworkflows:
        for unit_name in unit_names:
            if unit := subworkflow.get_unit_by_name(name=unit_name):
                patch_qe_input(unit, parameters, input_name=input_name)
                subworkflow.set_unit(unit)
