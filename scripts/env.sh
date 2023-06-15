#!/bin/bash

# To run as:
#   curl https://raw.githubusercontent.com/Exabyte-io/api-examples/${GIT_BRANCH}/scripts/env.sh | bash

set -euo pipefail

GIT_BRANCH="${GIT_BRANCH:-dev}"
REPO_NAME="api-examples"
VERBOSE="${EXABYTE_API_EXAMPLES_VERBOSE:-}"
IS_COLAB="${COLAB_JUPYTER_IP:-}"
NEED_GIT_LFS="${NEED_GIT_LFS:-}"

if [ ! -z "${VERBOSE}" ]; then
    set -vxeuo pipefail
    stdout="/dev/stdout"
else
    stdout="/dev/null"
fi

if [ ! -z "${IS_COLAB}" ]; then
    git clone https://github.com/Exabyte-io/${REPO_NAME}.git --single-branch --branch ${GIT_BRANCH} > ${stdout} 2>&1 || \
        echo -e "Directory ${REPO_NAME} already exists. Nothing to do." > ${stdout} 2>&1

    cd ${REPO_NAME}

    python -m pip install -r requirements-colab.txt > ${stdout} 2>&1
    python -m pip install examples/ > ${stdout} 2>&1

    if [ ! -z "${NEED_GIT_LFS}" ]; then
        sudo apt-get install -y git-lfs > ${stdout} 2>&1
        git lfs pull > ${stdout} 2>&1
    fi

    notebook_path="$(notebook-info)"  # comes from entry-points in 'setup.py'

    echo "Installation of the prerequisites is complete, the environment is ready to use!"
    echo -e "Notebook : ${notebook_path}"
fi
