import json
import os


async def authenticate_jupyterlite(data_from_host: dict) -> None:
    apiConfig = data_from_host.get("apiConfig")
    os.environ.update(data_from_host.get("environ", {}))
    os.environ.update(
        dict(
            ACCOUNT_ID=apiConfig.get("accountId"),  # type: ignore
            AUTH_TOKEN=apiConfig.get("authToken"),  # type: ignore
            ORGANIZATION_ID=apiConfig.get("organizationId", ""),  # type: ignore
            CLUSTERS=json.dumps(apiConfig.get("clusters", [])),  # type: ignore
        )
    )
