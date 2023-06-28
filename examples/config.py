import io
import os
import re

from notebook.utils import to_api_path

_script_exporter = None


def script_post_save(model, os_path, contents_manager, **kwargs):
    """
    Converts notebooks to Python script on save.
    See https://jupyter-notebook.readthedocs.io/en/stable/extending/savehooks.html for more information.
    """
    from nbconvert.exporters.script import ScriptExporter

    if model["type"] != "notebook":
        return

    global _script_exporter

    if _script_exporter is None:
        _script_exporter = ScriptExporter(parent=contents_manager)

    log = contents_manager.log

    base, ext = os.path.splitext(os.path.realpath(os_path))
    script, resources = _script_exporter.from_filename(os_path)
    script_fname = base + resources.get("output_extension", ".txt")
    log.info("Saving script /%s", to_api_path(script_fname, contents_manager.root_dir))

    # Remove notebook numberings from auto-generated .py files
    pattern = """
    (?<=^\#\sIn\[)   # Lookbehind to ensure we get the "# In[" bits notebooks put at the start of these lines
    (\d+)           # 1 or more digits, since we don't care about changing "[]" portions
    (?=\]:)         # Lookahead to get to the closing "]:" at the end of the line
    """
    script = re.sub(pattern, "", script, flags=re.MULTILINE | re.VERBOSE)

    with io.open(script_fname, "w", encoding="utf-8") as f:
        f.write(script)


c.FileContentsManager.post_save_hook = script_post_save  # noqa F821
