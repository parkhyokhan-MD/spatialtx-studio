from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree


def _upper_values(matrix: sparse.spmatrix) -> np.ndarray:
    upper = sparse.triu(matrix, k=1).tocoo()
    return np.asarray(upper.data, dtype=float)


def _coordinate_status(coords: np.ndarray | None, n_nodes: int) -> tuple[bool, list[str], int]:
    warnings: list[str] = []
    duplicate_count = 0
    if coords is None:
        return False, ["missing spatial coordinates"], 0
    coords = np.asarray(coords, dtype=float)
    valid = coords.ndim == 2 and coords.shape == (n_nodes, 2) and np.isfinite(coords).all()
    if not valid:
        warnings.append("non-finite or wrong-shape spatial coordinates")
        return False, warnings, 0
    if n_nodes:
        unique = np.unique(coords, axis=0)
        duplicate_count = int(n_nodes - len(unique))
        if duplicate_count:
            warnings.append(f"duplicate coordinates detected: {duplicate_count}")
    return True, warnings, duplicate_count


def graph_qc(
    connectivities: sparse.spmatrix,
    distances: sparse.spmatrix,
    coords: np.ndarray | None = None,
    *,
    method: str = "radius",
    typical_spacing: float | None = None,
    long_edge_factor: float = 3.0,
    isolated_fraction_warning: float = 0.10,
    largest_component_ratio_warning: float = 0.80,
    median_degree_warning: float = 2.0,
    long_edge_fraction_warning: float = 0.05,
    near_complete_density_warning: float = 0.50,
) -> tuple[dict, list[str], pd.DataFrame]:
    """Calculate graph QC metrics without densifying NxN matrices."""
    conn = connectivities.tocsr()
    dist = distances.tocsr()
    n_nodes = int(conn.shape[0])
    if conn.shape != dist.shape:
        raise ValueError("connectivities and distances must have the same shape")
    if conn.shape[0] != conn.shape[1]:
        raise ValueError("graph matrices must be square")
    undirected_edges = sparse.triu(conn, k=1).tocoo()
    n_edges = int(undirected_edges.nnz)
    degree = np.asarray((conn > 0).sum(axis=1)).ravel().astype(float)
    edge_distances = _upper_values(dist)
    n_components = int(connected_components(conn > 0, directed=False, return_labels=False)) if n_nodes else 0
    if n_nodes:
        component_count, component_labels = connected_components(conn > 0, directed=False, return_labels=True)
        sizes = np.bincount(component_labels, minlength=component_count) if component_count else np.asarray([], dtype=int)
        largest_component_ratio = float(sizes.max() / n_nodes) if len(sizes) else 0.0
    else:
        component_count = 0
        largest_component_ratio = 0.0
    coordinate_valid, coordinate_warnings, duplicate_count = _coordinate_status(coords, n_nodes)
    warnings = list(coordinate_warnings)
    isolated_fraction = float(np.mean(degree == 0)) if n_nodes else 0.0
    density = float((2 * n_edges) / (n_nodes * (n_nodes - 1))) if n_nodes > 1 else 0.0
    finite_dist = edge_distances[np.isfinite(edge_distances) & (edge_distances > 0)]
    if typical_spacing is None:
        typical_spacing = float(np.median(finite_dist)) if len(finite_dist) else np.nan
    long_threshold = float(typical_spacing * long_edge_factor) if np.isfinite(typical_spacing) else np.nan
    long_edge_fraction = (
        float(np.mean(finite_dist > long_threshold))
        if len(finite_dist) and np.isfinite(long_threshold)
        else 0.0
    )
    nearest_distances = np.asarray([], dtype=float)
    if coordinate_valid and n_nodes > 1:
        nearest_distances = np.asarray(cKDTree(np.asarray(coords, dtype=float)).query(coords, k=2)[0][:, 1], dtype=float)
        nearest_distances = nearest_distances[np.isfinite(nearest_distances)]
    if n_edges == 0:
        warnings.append("graph has no usable spatial edges")
    if isolated_fraction > float(isolated_fraction_warning):
        warnings.append(f"excessive isolated spots: {isolated_fraction:.1%}")
    if n_nodes and largest_component_ratio < float(largest_component_ratio_warning):
        warnings.append(f"largest connected component is too small: {largest_component_ratio:.1%}")
    median_degree = float(np.median(degree)) if n_nodes else 0.0
    if n_nodes and median_degree < float(median_degree_warning):
        warnings.append(f"low median degree: {median_degree:.3g}")
    if method == "knn" and long_edge_fraction > float(long_edge_fraction_warning):
        warnings.append(f"KNN created physically long edges: {long_edge_fraction:.1%}")
    if method == "radius" and n_nodes > 1:
        # Density remains descriptive only. A low density is expected for a
        # well-behaved local graph as n grows and is never a standalone warning.
        if density > float(near_complete_density_warning):
            warnings.append("radius graph is nearly complete")
    quantiles = [0.05, 0.25, 0.50, 0.75, 0.95]
    if len(finite_dist):
        q_values = np.quantile(finite_dist, quantiles)
    else:
        q_values = np.full(len(quantiles), np.nan)
    metrics = {
        "n_nodes": n_nodes,
        "n_edges": n_edges,
        "mean_degree": float(np.mean(degree)) if n_nodes else 0.0,
        "median_degree": median_degree,
        "min_degree": float(np.min(degree)) if n_nodes else 0.0,
        "max_degree": float(np.max(degree)) if n_nodes else 0.0,
        "isolated_fraction": isolated_fraction,
        "n_connected_components": int(n_components),
        "largest_component_ratio": largest_component_ratio,
        "edge_distance_median": float(np.median(finite_dist)) if len(finite_dist) else np.nan,
        "edge_distance_q05": float(q_values[0]),
        "edge_distance_q25": float(q_values[1]),
        "edge_distance_q50": float(q_values[2]),
        "edge_distance_q75": float(q_values[3]),
        "edge_distance_q95": float(q_values[4]),
        "long_edge_fraction": long_edge_fraction,
        "graph_density": density,
        "coordinate_valid": bool(coordinate_valid),
        "duplicate_coordinate_count": duplicate_count,
        "nearest_neighbor_distance_min": float(np.min(nearest_distances)) if len(nearest_distances) else np.nan,
        "nearest_neighbor_distance_median": float(np.median(nearest_distances)) if len(nearest_distances) else np.nan,
        "nearest_neighbor_distance_q95": float(np.quantile(nearest_distances, 0.95)) if len(nearest_distances) else np.nan,
        "number_of_spots": n_nodes,
        "number_of_edges": n_edges,
        "minimum_degree": float(np.min(degree)) if n_nodes else 0.0,
        "maximum_degree": float(np.max(degree)) if n_nodes else 0.0,
        "isolated_spot_fraction": isolated_fraction,
        "number_of_connected_components": int(n_components),
        "qc_threshold_isolated_spot_fraction": float(isolated_fraction_warning),
        "qc_threshold_largest_component_ratio": float(largest_component_ratio_warning),
        "qc_threshold_median_degree": float(median_degree_warning),
        "qc_threshold_long_edge_fraction": float(long_edge_fraction_warning),
        "qc_threshold_near_complete_density": float(near_complete_density_warning),
        "density_role": "informational_not_standalone_failure",
    }
    degree_table = pd.DataFrame({"spot_index": np.arange(n_nodes), "degree": degree.astype(int)})
    return metrics, warnings, degree_table
