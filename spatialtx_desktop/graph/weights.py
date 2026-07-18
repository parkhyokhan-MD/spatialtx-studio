from __future__ import annotations

import numpy as np


def estimate_radius(coords: np.ndarray, multiplier: float = 1.25) -> float:
    """Estimate a conservative radius in the active coordinate unit."""
    from scipy.spatial import cKDTree

    coords = np.asarray(coords, dtype=float)
    if len(coords) < 2:
        return 0.0
    distances, _ = cKDTree(coords).query(coords, k=min(2, len(coords)))
    nearest = np.asarray(distances[:, 1], dtype=float)
    nearest = nearest[np.isfinite(nearest) & (nearest > 0)]
    if not len(nearest):
        return 0.0
    return float(np.median(nearest) * multiplier)


def edge_weights(distances: np.ndarray, method: str = "binary", gaussian_sigma: float | None = None) -> np.ndarray:
    """Convert edge distances into connectivities."""
    distances = np.asarray(distances, dtype=float)
    method = (method or "binary").lower()
    if method == "binary":
        return np.ones_like(distances, dtype=float)
    if method == "inverse_distance":
        return 1.0 / np.maximum(distances, np.finfo(float).eps)
    if method == "gaussian":
        finite = distances[np.isfinite(distances) & (distances > 0)]
        sigma = float(gaussian_sigma or 0.0)
        if sigma <= 0:
            sigma = float(np.median(finite)) if len(finite) else 1.0
        sigma = max(sigma, np.finfo(float).eps)
        return np.exp(-0.5 * (distances / sigma) ** 2)
    raise ValueError("edge weighting must be one of: binary, inverse_distance, gaussian")
