from mat3ra.api_client import APIClient


def get_or_create_material(api_client: APIClient, material, owner_id: str) -> dict:
    """
    Returns an existing material from the collection if one with the same structural hash
    exists under the given owner, otherwise creates a new one.

    Args:
        api_client (APIClient): API client instance carrying the authorization context.
        material: mat3ra-made Material object (must have a .hash property).
        owner_id (str): Account ID under which to search and create.

    Returns:
        dict: The material dict (existing or newly created).
    """
    existing = api_client.materials.list({"hash": material.hash, "owner._id": owner_id})
    if existing:
        print(f"♻️  Reusing already existing Material: {existing[0]['_id']}")
        return existing[0]
    created = api_client.materials.create(material.to_dict(), owner_id=owner_id)
    print(f"✅ Material created: {created['_id']}")
    return created
