import asyncio
import os
import time

import requests
from IPython.display import Javascript, display

OIDC_BASE_URL = "http://localhost:3000/oidc"  # Your OIDC server URL
CLIENT_ID = "default-client"  # Your OAuth client ID
CLIENT_SECRET = "default-secret"  # Your OAuth client secret
SCOPE = "openid profile email"  # Requested scopes


async def authenticate_device_flow(
    oidc_base_url=OIDC_BASE_URL,
    client_id=CLIENT_ID,
    client_secret=CLIENT_SECRET,
    scope=SCOPE,
):
    device_response = requests.post(
        f"{oidc_base_url}/device/auth",
        data={"client_id": client_id, "client_secret": client_secret, "scope": scope},
        headers={"Content-Type": "application/x-www-form-urlencoded"},
        timeout=10,
    )
    device_response.raise_for_status()

    device_data = device_response.json()
    device_code = device_data["device_code"]
    user_code = device_data["user_code"]
    verification_uri_complete = device_data.get("verification_uri_complete", device_data["verification_uri"])
    polling_interval_seconds = int(device_data.get("interval", 5))
    expires_in_seconds = int(device_data.get("expires_in", 600))

    # JupyterLite: window.open must happen during cell execution to avoid popup blocker.
    display(
        Javascript(
            "alert("
            + repr(f"Open the login page and enter this code:\n\n{user_code}")
            + ");"
            + f"window.open({verification_uri_complete!r}, '_blank');"
        )
    )

    start_time_seconds = time.time()
    while time.time() - start_time_seconds < expires_in_seconds:
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
            token_data = token_response.json()
            os.environ["OIDC_ACCESS_TOKEN"] = token_data["access_token"]
            if "refresh_token" in token_data:
                os.environ["OIDC_REFRESH_TOKEN"] = token_data["refresh_token"]
            return token_data

        error_data = (
            token_response.json()
            if token_response.headers.get("content-type", "").startswith("application/json")
            else {}
        )
        error_code = error_data.get("error", "")

        if error_code == "slow_down":
            polling_interval_seconds += 5
        elif error_code != "authorization_pending":
            raise Exception(error_data.get("error_description") or error_code or token_response.text)

        await asyncio.sleep(polling_interval_seconds)

    raise Exception("Timeout waiting for authorization.")
