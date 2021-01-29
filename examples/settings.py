# 
# This file contains common variables that are used inside the examples.
#

# General settings. For how the notebooks operate.

    # use_interactive_JSON_viewer: Whether to use the IPython interactive viewer, or print in plaintext.

use_interactive_JSON_viewer = False

# Account settings. Need a one-time adjustment for examples to work.

    # ACCOUNT_ID: Account ID. See get_authentication_params.ipynb example for more information.
    # AUTH_TOKEN: Account authentication token. See get_authentication_params.ipynb for more information.
    # MATERIALS_PROJECT_API_KEY: Materials project API key. See https://materialsproject.org/open for more information.

ACCOUNT_ID = "ACCOUNT_ID"
AUTH_TOKEN = "AUTH_TOKEN"
MATERIALS_PROJECT_API_KEY = "MATERIALS_PROJECT_API_KEY"

# Advanced settings. Should not need adjustments.

    # HOST: Hostname of the RESTful API server. Defaults to platform.exabyte.io.
    # PORT: The port RESTful API server is listening on. Defaults to 443.
    # VERSION: RESTFul API version. Defaults to 2018-10-01.
    # SECURE: Whether to use secure connection. Defaults to True.

PORT = 443
SECURE = True
VERSION = "2018-10-01"
HOST = "platform.exabyte.io"
ENDPOINT_ARGS = [HOST, PORT, ACCOUNT_ID, AUTH_TOKEN, VERSION, SECURE]
