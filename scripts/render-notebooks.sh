#!/bin/bash

# This is a helper script to render Jupyter notebook files (.ipynb)
#
# The workflow is as follows:
# - parse the README.md file and extract all the notebooks listed in the table
#   or use a list of notebooks passed as a blank-separated inputs in command line;
# - for each notebook:
#   - execute `nbstripout` to strip output cells;
#   - execute `nbconvert` to populate the .ipynb file in place;
#   - execute `nbconvert` to .py.

set -euo pipefail

if [ ! -z "$*" ]; then
    notebooks="$@"
    # E.g.:
    # notebooks="examples/material/api_interoperability_showcase.ipynb"
else
    notebooks="$(cat README.md | grep -E "\| (.*)[(.*)](examples/(.*).ipynb)" | cut -d'|' -f3 | grep -oE "\((.*)\)" | cut -d'(' -f2 | cut -d')' -f1)"
fi

idir="$PWD"

function now() {
     echo $(date '+%F %T')
}

for notebook in ${notebooks}; do
    echo -e "$(now) Processing ${notebook}..."

    if [ ! -f "${notebook}" ]; then
        echo "File ${notebook} does not exist. Exiting."
        exit 1
    fi

    notebook_dir=$(dirname $notebook)
    notebook_name=$(basename $notebook)

    cd $notebook_dir

    echo -e "$(now) Stripout output cells (if any)"
    nbstripout ${notebook_name}

    # Convert stripout notebooks to Python
    echo -e "$(now) Exporting a python script for ${notebook_name}..."
    jupyter-nbconvert --to python $notebook_name
    echo -e "$(now) Exporting of a python script for ${notebook_name} is complete."

    # Execute and save notebooks to HTML
    if [[ ("${notebook_name}" == *"get_authentication_params.ipynb") || \
          ("${notebook_name}" == *"run-simulations-and-extract-properties.ipynb") || \
          ("${notebook_name}" == *"this-notebook-does-not-exist--placeholder-for-future-updates.ipynb") \
          ]]; then
        echo -e "$(now) Exporting to html without execution of ${notebook_name}..."
        jupyter-nbconvert --to html ${notebook_name}
        echo -e "$(now) Exporting to html without execution of ${notebook_name} is complete."
    else
        echo -e "$(now) Executing ${notebook_name} and saving to html format..."
        jupyter-nbconvert --execute --to html ${notebook_name}
        echo -e "$(now) Execution of ${notebook_name} is complete."
    fi

    cd $idir
    echo ""
done

exit 0
