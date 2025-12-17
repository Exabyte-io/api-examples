import asyncio
import os
import time

import requests
from IPython.display import Javascript, display

# Default OIDC Configuration
OIDC_BASE_URL = "http://localhost:3000/oidc"  # Your OIDC server URL
CLIENT_ID = "default-client"  # Your OAuth client ID
CLIENT_SECRET = "default-secret"  # Your OAuth client secret
SCOPE = "openid profile email"  # Requested scopes

# Environment Variable Names
ACCESS_TOKEN_ENV_VAR = "OIDC_ACCESS_TOKEN"
REFRESH_TOKEN_ENV_VAR = "OIDC_REFRESH_TOKEN"


def request_device_flow_state(oidc_base_url: str, client_id: str, client_secret: str, scope: str) -> dict:
    device_response = requests.post(
        f"{oidc_base_url}/device/auth",
        data={"client_id": client_id, "client_secret": client_secret, "scope": scope},
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
    display(
        Javascript(
            "alert("
            + repr(f"Open the login page and enter this code:\\n\\n{user_code}")
            + ");"
            + f"window.open({verification_uri_complete!r}, '_blank');"
        )
    )


def store_token_data_in_environment(token_data: dict) -> None:
    os.environ[ACCESS_TOKEN_ENV_VAR] = token_data["access_token"]
    if "refresh_token" in token_data:
        os.environ[REFRESH_TOKEN_ENV_VAR] = token_data["refresh_token"]


async def _poll_for_token_data(
    oidc_base_url: str,
    client_id: str,
    client_secret: str,
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
                "client_secret": client_secret,
            },
            headers={"Content-Type": "application/x-www-form-urlencoded"},
            timeout=10,
        )
        if token_response.status_code == 200:
            return token_response.json()
        await asyncio.sleep(polling_interval_seconds)
    raise Exception("Timeout waiting for authorization.")


async def authenticate(
    oidc_base_url=OIDC_BASE_URL,
    client_id=CLIENT_ID,
    client_secret=CLIENT_SECRET,
    scope=SCOPE,
) -> dict:
    device_flow_state = request_device_flow_state(oidc_base_url, client_id, client_secret, scope)
    show_device_flow_popup(device_flow_state["verification_uri_complete"], device_flow_state["user_code"])
    token_data = await _poll_for_token_data(
        oidc_base_url=oidc_base_url,
        client_id=client_id,
        client_secret=client_secret,
        device_code=device_flow_state["device_code"],
        polling_interval_seconds=device_flow_state["polling_interval_seconds"],
        expires_in_seconds=device_flow_state["expires_in_seconds"],
    )
    store_token_data_in_environment(token_data)
    return token_data
