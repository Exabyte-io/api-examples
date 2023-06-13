#!/bin/bash

# To run as:
#   curl https://raw.githubusercontent.com/Exabyte-io/api-examples/${GIT_BRANCH}/scripts/env.sh | bash

set -vxeuo pipefail

GIT_BRANCH="${GIT_BRANCH:-'dev'}"
REPO_NAME="api-examples"

git clone https://github.com/Exabyte-io/${REPO_NAME}.git --single-branch --branch ${GIT_BRANCH} || \
    echo -e "Directory ${REPO_NAME} already exists. Nothing to do."

python -m pip install -r ${REPO_NAME}/requirements-colab.txt
python -m pip install ./${REPO_NAME}/examples/
