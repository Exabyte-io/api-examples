# General settings. For how the notebooks operate.

from .python.other.api.settings import ACCOUNT_ID, AUTH_TOKEN, MATERIALS_PROJECT_API_KEY, ORGANIZATION_ID

# use_interactive_JSON_viewer: Whether to use the IPython interactive viewer, or print in plaintext.
use_interactive_JSON_viewer = True


UPLOADS_FOLDER = "uploads"

__all__ = ["ACCOUNT_ID", "AUTH_TOKEN", "MATERIALS_PROJECT_API_KEY", "ORGANIZATION_ID"]
