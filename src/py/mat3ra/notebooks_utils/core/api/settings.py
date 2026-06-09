#
# This file contains common variables that are used inside the examples.
#
import json
import os

# Account settings. Need a one-time adjustment for examples to work. These should be set in settings.json

# ACCOUNT_ID: Account ID. See get_authentication_params.ipynb example for more information.
# AUTH_TOKEN: Account authentication token. See get_authentication_params.ipynb for more information.
# MATERIALS_PROJECT_API_KEY: Materials project API key.
# See https://legacy.materialsproject.org/open for more information.

# Load variables from the settings.json file

absolute_path_to_settings_json_file = os.path.join(os.path.dirname(__file__), "settings.json")
assert absolute_path_to_settings_json_file
with open(absolute_path_to_settings_json_file) as settings_json_file:
    settings_json_config = json.load(settings_json_file)

# Environment variables take precedence over settings.json (set by authenticate() or JupyterLite host).
# For local Jupyter, set legacy API token values in settings.json or use OIDC via authenticate().

ACCOUNT_ID = os.getenv("ACCOUNT_ID", settings_json_config.get("ACCOUNT_ID"))
AUTH_TOKEN = os.getenv("AUTH_TOKEN", settings_json_config.get("AUTH_TOKEN"))
MATERIALS_PROJECT_API_KEY = os.getenv(
    "MATERIALS_PROJECT_API_KEY", settings_json_config.get("MATERIALS_PROJECT_API_KEY")
)
ORGANIZATION_ID = os.getenv("ORGANIZATION_ID", settings_json_config.get("ORGANIZATION_ID"))

# Advanced settings. Should not need adjustments.

# HOST: Hostname of the RESTful API server. Defaults to platform.mat3ra.com.
# PORT: The port RESTful API server is listening on. Defaults to 443.
# VERSION: RESTFul API version. Defaults to 2018-10-01.
# SECURE: Whether to use secure connection. Defaults to True.

PORT = os.environ.get("API_PORT", 443)
SECURE = os.environ.get("API_SECURE") != "false"
VERSION = os.environ.get("API_VERSION", "2018-10-01")
HOST = os.environ.get("API_HOST", "platform.mat3ra.com")
ENDPOINT_ARGS = [HOST, PORT, ACCOUNT_ID, AUTH_TOKEN, VERSION, SECURE]
