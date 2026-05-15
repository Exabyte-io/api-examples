import asyncio
import os
import time
from typing import Optional

import requests
from IPython.display import Javascript, display
from mat3ra.api_client import ACCESS_TOKEN_ENV_VAR, CLIENT_ID, SCOPE, APIEnv, build_oidc_base_url

REFRESH_TOKEN_ENV_VAR = "OIDC_REFRESH_TOKEN"


def get_oidc_base_url() -> str:
    env = APIEnv.from_env()
    return build_oidc_base_url(env.host, env.port, env.secure)


def request_device_flow_state(oidc_base_url: str, client_id: str, scope: str) -> dict:
    """
    Request an OAuth/OIDC Device Authorization flow state.

    Args:
        oidc_base_url: Base OIDC URL.
        client_id: OAuth client identifier for the device flow.
        scope: Space-separated scopes to request (e.g. "openid profile email").

    Returns:
        A dict with device_code, user_code, verification_uri_complete,
        polling_interval_seconds, expires_in_seconds.
    """
    device_response = requests.post(
        f"{oidc_base_url}/device/auth",
        data={"client_id": client_id, "scope": scope},
        headers={"Content-Type": "application/x-www-form-urlencoded"},
        timeout=10,
    )
    device_response.raise_for_status()
    device_data = device_response.json()
    return {
        "device_code": device_data["device_code"],
        "user_code": device_data["user_code"],
        "verification_uri_complete": device_data.get("verification_uri_complete", device_data["verification_uri"]),
        "polling_interval_seconds": int(device_data.get("interval", 5)),
        "expires_in_seconds": int(device_data.get("expires_in", 600)),
    }


def show_device_flow_popup(verification_uri_complete: str, user_code: str) -> None:
    from IPython.display import HTML

    display(
        HTML(
            f"<div style='padding: 15px; background: #e3f2fd; border-left: 4px solid #2196f3; margin: 10px 0;'>"
            f"<strong>Authentication Required</strong><br/>"
            f"Enter this code: <strong style='font-size: 1.2em; color: #1976d2;'>{user_code}</strong>"
            f"</div>"
        )
    )
    display(Javascript(f"window.open({verification_uri_complete!r}, '_blank');"))


def store_token_data_in_environment(token_data: dict) -> None:
    os.environ[ACCESS_TOKEN_ENV_VAR] = token_data["access_token"]
    if "refresh_token" in token_data:
        os.environ[REFRESH_TOKEN_ENV_VAR] = token_data["refresh_token"]


async def _poll_for_token_data(
    oidc_base_url: str,
    client_id: str,
    device_code: str,
    polling_interval_seconds: int,
    expires_in_seconds: int,
) -> dict:
    deadline_seconds = time.time() + expires_in_seconds
    while time.time() < deadline_seconds:
        token_response = requests.post(
            f"{oidc_base_url}/token",
            data={
                "grant_type": "urn:ietf:params:oauth:grant-type:device_code",
                "device_code": device_code,
                "client_id": client_id,
                "redirect_uris": [],
                "response_types": [],
                "token_endpoint_auth_method": "none",
            },
            headers={"Content-Type": "application/x-www-form-urlencoded"},
            timeout=10,
        )
        if token_response.status_code == 200:
            return token_response.json()
        await asyncio.sleep(polling_interval_seconds)
    raise Exception("Timeout waiting for authorization.")


async def authenticate_oidc(
    oidc_base_url: Optional[str] = None,
    client_id: str = CLIENT_ID,
    scope: str = SCOPE,
) -> dict:
    if oidc_base_url is None:
        oidc_base_url = get_oidc_base_url()
    device_flow_state = request_device_flow_state(oidc_base_url, client_id, scope)
    show_device_flow_popup(device_flow_state["verification_uri_complete"], device_flow_state["user_code"])
    token_data = await _poll_for_token_data(
        oidc_base_url=oidc_base_url,
        client_id=client_id,
        device_code=device_flow_state["device_code"],
        polling_interval_seconds=device_flow_state["polling_interval_seconds"],
        expires_in_seconds=device_flow_state["expires_in_seconds"],
    )
    store_token_data_in_environment(token_data)
    return token_data
