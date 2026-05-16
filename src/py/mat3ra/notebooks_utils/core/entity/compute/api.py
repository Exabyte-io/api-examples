import json
import os


def get_cluster_name(name: str = "cluster-001") -> str:
    clusters = json.loads(os.environ.get("CLUSTERS", "[]") or "[]")
    return clusters[0] if clusters else name
