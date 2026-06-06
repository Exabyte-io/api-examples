import io
import re
from typing import Any, Dict, List, Mapping, Optional, Tuple

import f90nml

JINJA_PATTERN = re.compile(r"\{\{.*?\}\}|\{%-?.*?-?%\}", re.DOTALL)


def _normalize_section(section: str) -> str:
    return section.lstrip("&").upper()


def _protect_jinja(text: str) -> Tuple[str, Dict[str, str]]:
    placeholders: Dict[str, str] = {}

    def replacer(match: re.Match[str]) -> str:
        key = f"__JINJA_{len(placeholders)}__"
        placeholders[key] = match.group(0)
        return f"'{key}'"

    return JINJA_PATTERN.sub(replacer, text), placeholders


def _restore_jinja(text: str, placeholders: Mapping[str, str]) -> str:
    restored = text
    for key, original in placeholders.items():
        restored = restored.replace(f"'{key}'", original)
    return restored


def _extract_namelist_block(content: str, section: str) -> Tuple[str, str, str]:
    normalized = _normalize_section(section)
    pattern = rf"(?ms)(^&{re.escape(normalized)}\s*\n.*?^/\s*$)"
    match = re.search(pattern, content)
    if match is None:
        raise ValueError(f"Namelist '&{normalized}' not found in input template.")
    return content[: match.start()], match.group(0), content[match.end() :]


def _restore_namelist_header(block: str, section: str) -> str:
    normalized = _normalize_section(section)
    return re.sub(r"^&\w+", f"&{normalized}", block, count=1, flags=re.MULTILINE)


def _patch_namelist_block(block: str, section: str, parameters: Mapping[str, Any]) -> str:
    normalized = _normalize_section(section)
    protected, placeholders = _protect_jinja(block)
    nml = f90nml.reads(protected)
    nml.patch({normalized.lower(): dict(parameters)})
    buffer = io.StringIO()
    nml.write(buffer)
    patched = _restore_namelist_header(buffer.getvalue(), normalized)
    return _restore_jinja(patched, placeholders)


def set_content(content: str, section: str, parameters: Mapping[str, Any]) -> str:
    before, block, after = _extract_namelist_block(content, section)
    patched_block = _patch_namelist_block(block, section, parameters)
    return before + patched_block + after


def _get_template(input_item):
    if isinstance(input_item, dict):
        template = input_item.get("template")
        if isinstance(template, dict):
            return template
        return input_item
    template = getattr(input_item, "template", None)
    if template is not None:
        return template
    return input_item


def _get_input_content(input_item) -> Optional[str]:
    template = _get_template(input_item)
    if isinstance(template, dict):
        return template.get("content")
    return getattr(template, "content", None)


def _set_input_content(input_item, content: str) -> None:
    template = _get_template(input_item)
    if isinstance(template, dict):
        template["content"] = content
    else:
        template.content = content


def _get_input_name(input_item) -> Optional[str]:
    template = _get_template(input_item)
    if isinstance(template, dict):
        return template.get("name")
    return getattr(template, "name", None)


def _is_multi_section_parameters(value: Any) -> bool:
    if not isinstance(value, Mapping):
        return False
    return all(isinstance(item, Mapping) for item in value.values())


def _patch_qe_input_single(
    unit,
    section: str,
    parameters: Mapping[str, Any],
    input_name: Optional[str] = None,
) -> None:
    matching_inputs = 0
    for input_item in getattr(unit, "input", []):
        content = _get_input_content(input_item)
        if content is None:
            continue
        if input_name is not None and _get_input_name(input_item) != input_name:
            continue
        _set_input_content(input_item, set_content(content, section, parameters))
        matching_inputs += 1

    if matching_inputs == 0:
        raise ValueError("No matching input template found for QE patch.")


def patch_qe_input(
    unit,
    section_or_parameters: Any,
    parameters: Optional[Mapping[str, Any]] = None,
    input_name: Optional[str] = None,
) -> None:
    if parameters is None and _is_multi_section_parameters(section_or_parameters):
        for section, section_parameters in section_or_parameters.items():
            _patch_qe_input_single(unit, section, section_parameters, input_name)
        return
    if parameters is None:
        raise TypeError("Expected section parameters mapping or a section name with parameters.")
    _patch_qe_input_single(unit, section_or_parameters, parameters, input_name)


def patch_workflow_qe_input(
    workflow,
    parameters: Mapping[str, Mapping[str, Any]],
    unit_names: List[str],
    input_name: Optional[str] = None,
) -> None:
    for subworkflow in workflow.subworkflows:
        for unit_name in unit_names:
            unit = subworkflow.get_unit_by_name(name=unit_name)
            if unit:
                patch_qe_input(unit, parameters, input_name=input_name)
                subworkflow.set_unit(unit)
