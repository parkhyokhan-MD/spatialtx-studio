from __future__ import annotations


def resolve_distance_config(config: dict) -> dict:
    distance = config.get("distance", {})
    mode = distance.get("mode", "spot")
    if mode != "spot":
        return {
            **distance,
            "mode": "spot",
            "note": "Physical distance calibration is reserved for a future version; v0.1 uses spot distance.",
        }
    return {**distance, "note": "Using spot-based distance."}
