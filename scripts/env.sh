#!/bin/bash

# To run as:
#   curl https://raw.githubusercontent.com/Exabyte-io/api-examples/${GIT_BRANCH}/scripts/env.sh | bash

set -vxeuo pipefail

GIT_BRANCH="${GIT_BRANCH:-'dev'}"

git clone https://github.com/Exabyte-io/api-examples.git --single-branch --branch ${GIT_BRANCH}

python -m pip install -r api-examples/requirements-colab.txt
python -m pip install ./api-examples/utils/

# python api-examples/examples/utils/initialize_settings.py  # the script with cd into the current notebook dir
