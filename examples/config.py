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

    if model['type'] != 'notebook':
        return

    global _script_exporter

    if _script_exporter is None:
        _script_exporter = ScriptExporter(parent=contents_manager)

    log = contents_manager.log

    base, ext = os.path.splitext(os.path.realpath(os_path))
    script, resources = _script_exporter.from_filename(os_path)
    script_fname = base + resources.get('output_extension', '.txt')
    log.info("Saving script /%s", to_api_path(script_fname, contents_manager.root_dir))

    # Because we now install modules in Notebooks, jupyter will save things with extra get_ipython() calls, that
    # regular python doesn't have. Replace those with normal calls to pip before saving the .py file.
    pattern = """
    ^                                   # Beginning of string
    (\s+)                               # Capture group 1, to preserve indentation
    .*                                  # The get_ipython() parts that we're removing
    (?<=get_ipython\(\).system)         # Lookbehind to ensure we only get the get_ipython lines
        \('                             # Beginning parenthesis and string quote
            \{sys.executable\}          # Call to sys.executable
            \s*-m\s*pip\s*install\s*    # Pip install command
            (.*)                        # Capture group 2, name of the module we want to install
        '\)                             # Ending parenthesis and string quote
    \s*                                 # Whitespace at the end of the line
    $                                   # End of string
    """
    # Replace get_ipython stuff with calls that work outside ipy notebooks. For example, this:
    #   get_ipython().system('{sys.executable} -m pip install pandas==1.1.4')
    # Becomes this:
    #   import subprocess, sys
    #   subprocess.call([sys.executable,'-m','pip','install','pandas==1.1.4'])
    script = re.sub(pattern, r"\1import subprocess, sys\n\1subprocess.call([sys.executable,'-m','pip','install','\2'])",
                    script, flags=re.MULTILINE | re.VERBOSE)

    with io.open(script_fname, 'w', encoding='utf-8') as f:
        f.write(script)


c.FileContentsManager.post_save_hook = script_post_save
