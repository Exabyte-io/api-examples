#!/bin/bash

# This is a helper script to render Jupyter notebook files (.ipynb)
#
# The workflow is as follows:
# - parse the README.md file and extract all the notebooks listed in the table;
# - for each notebook:
#   - execute `nbstripout` to strip output cells;
#   - execute `nbconvert` to populate the .ipynb file in place;
#   - execute `nbconvert` to .py.

set -euo pipefail

notebooks="$(cat README.md | grep -E "\| (.*)[(.*)](examples/(.*).ipynb)" | cut -d'|' -f3 | grep -oE "\((.*)\)" | cut -d'(' -f2 | cut -d')' -f1)"

idir="$PWD"

function now() {
     echo $(date '+%F %T')
}

for notebook in ${notebooks}; do
    echo -e "$(now) Processing ${notebook}..."

    notebook_dir=$(dirname $notebook)
    notebook_name=$(basename $notebook)

    cd $notebook_dir

    echo -e "$(now) Stripout output cells (if any)"
    nbstripout ${notebook_name}

    if [[ ("${notebook_name}" == *"get_authentication_params.ipynb") || \
          ("${notebook_name}" == *"run-simulations-and-extract-properties.ipynb") || \
          ("${notebook_name}" == *"this-notebook-does-not-exist--placeholder-for-future-updates.ipynb") \
          ]]; then
        echo -e "$(now) Skipping execution of ${notebook_name}."
    else
        echo -e "$(now) Executing ${notebook_name} in place..."
        jupyter-nbconvert --execute --inplace ${notebook_name}
        echo -e "$(now) Execution of ${notebook_name} in complete."
    fi

    echo -e "$(now) Exporting a python script for ${notebook_name}..."
    jupyter-nbconvert --to python $notebook_name
    echo -e "$(now) Exporting of a python script for ${notebook_name} is complete."

    cd $idir
    echo ""
done
