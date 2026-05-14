from typing import Any

from mat3ra.made.material import Material


def get_bulk_material(api_client: Any, slab_material: Material, owner_id: str):
    slab_dict = slab_material.to_dict()
    metadata = slab_dict.get("metadata") or {}
    bulk_crystal = None

    if metadata.get("bulkId") is not None:
        bulk_query = {"_id": metadata["bulkId"]}
    else:
        for build_step in reversed(metadata.get("build") or []):
            try:
                bulk_crystal = build_step["configuration"]["stack_components"][0]["crystal"]
                break
            except (KeyError, IndexError, TypeError):
                continue

        if bulk_crystal is None:
            raise ValueError(
                "No metadata.build[*].configuration.stack_components[0].crystal entry was found on the slab."
            )

        if bulk_crystal.get("_id") is not None:
            bulk_query = {"_id": bulk_crystal["_id"]}
        elif bulk_crystal.get("scaledHash") is not None:
            bulk_query = {"scaledHash": bulk_crystal["scaledHash"]}
        elif bulk_crystal.get("hash") is not None:
            bulk_query = {"hash": bulk_crystal["hash"]}
        else:
            try:
                bulk_query = {"hash": Material.create(bulk_crystal).hash}
            except Exception as exc:
                raise ValueError("Could not resolve a bulk query from the slab metadata.") from exc

    matches = api_client.materials.list(bulk_query)
    bulk_material_response = next(
        (item for item in matches if item.get("owner", {}).get("_id") == owner_id),
        None,
    ) or (matches[0] if matches else None)

    if bulk_material_response is None:
        raise ValueError(
            "The bulk material resolved from slab metadata is not present on the platform. "
            "Run the Total Energy notebook for that bulk material first, then rerun this notebook."
        )

    print(f"Found exact bulk material: {bulk_material_response['_id']}")

    return bulk_query, bulk_material_response, Material.create(bulk_material_response)
