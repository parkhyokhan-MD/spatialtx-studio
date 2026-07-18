from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Any

import numpy as np
from scipy import sparse
from scipy.spatial import cKDTree

from .metadata import GRAPH_KEYS, GRAPH_SCHEMA_VERSION, json_safe
from .qc import graph_qc
from .weights import edge_weights, estimate_radius


@dataclass
class GraphBuildConfig:
    method: str = "radius"
    radius: float | None = None
    k: int = 6
    lattice_ring: int = 1
    weighting: str = "binary"
    gaussian_sigma: float | None = None
    symmetrization: str = "union"
    coordinate_source: str = "obsm/spatial"
    coordinate_unit: str = "native"
    source_coordinate_unit: str = ""
    coordinate_scale: float | None = None
    scale_source: str = ""
    platform_guess: str = "unknown"
    random_seed: int = 0
    long_edge_factor: float = 3.0
    isolated_fraction_warning: float = 0.10
    largest_component_ratio_warning: float = 0.80
    median_degree_warning: float = 2.0
    long_edge_fraction_warning: float = 0.05
    near_complete_density_warning: float = 0.50


@dataclass
class GraphBuildResult:
    method: str
    connectivities: sparse.csr_matrix
    distances: sparse.csr_matrix
    metadata: dict[str, Any]
    qc: dict[str, Any]
    warnings: list[str] = field(default_factory=list)
    degree_table: Any | None = None


def _extract_spatial(adata) -> tuple[np.ndarray, str, list[str]]:
    warnings: list[str] = []
    if not hasattr(adata, "obsm") or "spatial" not in adata.obsm:
        return np.empty((0, 2), dtype=float), "obsm/spatial", ["missing spatial coordinates"]
    try:
        coords = np.asarray(adata.obsm["spatial"], dtype=float)
    except Exception:
        return np.empty((0, 2), dtype=float), "obsm/spatial", ["spatial coordinates are not numeric"]
    n_obs = int(getattr(adata, "n_obs", coords.shape[0] if coords.ndim else 0))
    if coords.ndim != 2 or coords.shape != (n_obs, 2):
        return coords, "obsm/spatial", ["spatial coordinates have the wrong shape"]
    if not np.isfinite(coords).all():
        warnings.append("spatial coordinates contain non-finite values")
    return coords, "obsm/spatial", warnings


def _empty_graph(n: int) -> tuple[sparse.csr_matrix, sparse.csr_matrix]:
    return sparse.csr_matrix((n, n), dtype=float), sparse.csr_matrix((n, n), dtype=float)


def _to_symmetric_sparse(
    n: int,
    rows: list[int],
    cols: list[int],
    distances: list[float],
    weighting: str,
    sigma: float | None,
) -> tuple[sparse.csr_matrix, sparse.csr_matrix, int]:
    if not rows:
        return *_empty_graph(n), 0
    excluded_zero_distance = 0
    if str(weighting).lower() == "inverse_distance":
        keep = np.asarray(distances, dtype=float) > 0
        excluded_zero_distance = int((~keep).sum())
        rows = [value for value, include in zip(rows, keep) if include]
        cols = [value for value, include in zip(cols, keep) if include]
        distances = [value for value, include in zip(distances, keep) if include]
        if not rows:
            return *_empty_graph(n), excluded_zero_distance
    row_array = np.asarray(rows + cols, dtype=int)
    col_array = np.asarray(cols + rows, dtype=int)
    dist_array = np.asarray(distances + distances, dtype=float)
    weights = edge_weights(dist_array, weighting, sigma)
    connectivities = sparse.coo_matrix((weights, (row_array, col_array)), shape=(n, n)).tocsr()
    distance_matrix = sparse.coo_matrix((dist_array, (row_array, col_array)), shape=(n, n)).tocsr()
    connectivities.eliminate_zeros()
    distance_matrix.eliminate_zeros()
    return connectivities, distance_matrix, excluded_zero_distance


def _radius_graph(coords: np.ndarray, radius: float | None, config: GraphBuildConfig) -> tuple[sparse.csr_matrix, sparse.csr_matrix, float, list[str], int]:
    warnings: list[str] = []
    n = len(coords)
    if n < 2:
        return *_empty_graph(n), 0.0, ["radius graph requires at least two spots"], 0
    radius_value = float(radius) if radius is not None and float(radius) > 0 else estimate_radius(coords)
    if radius_value <= 0:
        return *_empty_graph(n), radius_value, ["unable to estimate a positive graph radius"], 0
    tree = cKDTree(coords)
    pairs = sorted(tree.query_pairs(radius_value))
    rows: list[int] = []
    cols: list[int] = []
    distances: list[float] = []
    for i, j in pairs:
        rows.append(int(i))
        cols.append(int(j))
        distances.append(float(np.linalg.norm(coords[i] - coords[j])))
    conn, dist, excluded = _to_symmetric_sparse(n, rows, cols, distances, config.weighting, config.gaussian_sigma)
    return conn, dist, radius_value, warnings, excluded


def _knn_graph(coords: np.ndarray, config: GraphBuildConfig) -> tuple[sparse.csr_matrix, sparse.csr_matrix, list[str], int]:
    warnings: list[str] = []
    n = len(coords)
    if n < 2:
        return *_empty_graph(n), ["KNN graph requires at least two spots"], 0
    k = max(1, min(int(config.k), n - 1))
    tree = cKDTree(coords)
    distances, indices = tree.query(coords, k=k + 1)
    directed: set[tuple[int, int]] = set()
    distance_lookup: dict[tuple[int, int], float] = {}
    for i in range(n):
        for distance, j in zip(np.ravel(distances[i])[1:], np.ravel(indices[i])[1:]):
            if int(j) == i or not np.isfinite(distance):
                continue
            directed.add((i, int(j)))
            distance_lookup[(i, int(j))] = float(distance)
    rows: list[int] = []
    cols: list[int] = []
    edge_distances: list[float] = []
    seen: set[tuple[int, int]] = set()
    for i, j in sorted(directed):
        a, b = sorted((i, j))
        if (a, b) in seen:
            continue
        if config.symmetrization == "mutual" and not ((i, j) in directed and (j, i) in directed):
            continue
        seen.add((a, b))
        d1 = distance_lookup.get((a, b), np.inf)
        d2 = distance_lookup.get((b, a), np.inf)
        rows.append(a)
        cols.append(b)
        edge_distances.append(float(min(d1, d2)))
    conn, dist, excluded = _to_symmetric_sparse(n, rows, cols, edge_distances, config.weighting, config.gaussian_sigma)
    return conn, dist, warnings, excluded


def _lattice_offsets(ring: int) -> set[tuple[int, int]]:
    immediate = {(-1, -1), (-1, 1), (0, -2), (0, 2), (1, -1), (1, 1)}
    ring = max(1, int(ring))
    if ring == 1:
        return immediate
    offsets: set[tuple[int, int]] = set()
    frontier = {(0, 0)}
    visited = {(0, 0)}
    for _ in range(ring):
        next_frontier: set[tuple[int, int]] = set()
        for base_r, base_c in frontier:
            for dr, dc in immediate:
                item = (base_r + dr, base_c + dc)
                if item not in visited:
                    visited.add(item)
                    next_frontier.add(item)
                    offsets.add(item)
        frontier = next_frontier
    offsets.discard((0, 0))
    return offsets


def _as_obs_array(adata, name: str) -> np.ndarray | None:
    if not hasattr(adata, "obs") or name not in adata.obs:
        return None
    values = np.asarray(adata.obs[name])
    try:
        values = values.astype(float)
    except Exception:
        return None
    if not np.isfinite(values).all():
        return None
    return values


def _lattice_graph(
    adata,
    coords: np.ndarray,
    config: GraphBuildConfig,
) -> tuple[sparse.csr_matrix, sparse.csr_matrix, bool, str, float | None, list[str], int]:
    warnings: list[str] = []
    n = len(coords)
    rows_arr = _as_obs_array(adata, "array_row")
    cols_arr = _as_obs_array(adata, "array_col")
    if rows_arr is None or cols_arr is None or len(rows_arr) != n or len(cols_arr) != n:
        warnings.append("Visium lattice metadata unavailable; falling back to radius graph")
        conn, dist, radius_value, radius_warnings, excluded = _radius_graph(coords, config.radius, config)
        return conn, dist, False, "radius", radius_value, warnings + radius_warnings, excluded
    if not (np.allclose(rows_arr, np.round(rows_arr)) and np.allclose(cols_arr, np.round(cols_arr))):
        warnings.append("Visium lattice metadata is not integer-like; falling back to radius graph")
        conn, dist, radius_value, radius_warnings, excluded = _radius_graph(coords, config.radius, config)
        return conn, dist, False, "radius", radius_value, warnings + radius_warnings, excluded
    positions = {(int(r), int(c)): i for i, (r, c) in enumerate(zip(rows_arr, cols_arr))}
    if len(positions) != n:
        warnings.append("duplicate Visium lattice positions detected")
    offsets = _lattice_offsets(config.lattice_ring)
    rows: list[int] = []
    cols: list[int] = []
    distances: list[float] = []
    seen: set[tuple[int, int]] = set()
    for i, (r, c) in enumerate(zip(rows_arr.astype(int), cols_arr.astype(int))):
        for dr, dc in offsets:
            j = positions.get((int(r + dr), int(c + dc)))
            if j is None or j == i:
                continue
            a, b = sorted((i, int(j)))
            if (a, b) in seen:
                continue
            seen.add((a, b))
            rows.append(a)
            cols.append(b)
            distances.append(float(np.linalg.norm(coords[a] - coords[b])) if len(coords) == n else 1.0)
    if not rows:
        warnings.append("Visium lattice detection succeeded but produced no edges; falling back to radius graph")
        conn, dist, radius_value, radius_warnings, excluded = _radius_graph(coords, config.radius, config)
        return conn, dist, False, "radius", radius_value, warnings + radius_warnings, excluded
    conn, dist, excluded = _to_symmetric_sparse(n, rows, cols, distances, config.weighting, config.gaussian_sigma)
    return conn, dist, True, "lattice", None, warnings, excluded


def build_spatial_graph(adata, config: GraphBuildConfig | None = None) -> GraphBuildResult:
    """Build a Sparse SpatialTX graph from AnnData coordinates."""
    config = config or GraphBuildConfig()
    requested_method = (config.method or "radius").lower()
    if requested_method not in {"radius", "lattice", "knn"}:
        raise ValueError("graph method must be one of: radius, lattice, knn")
    coords, coord_source, coord_warnings = _extract_spatial(adata)
    audit = dict(getattr(adata, "uns", {}).get("spatialtx_input_audit", {}))
    config = replace(
        config,
        platform_guess=(str(audit.get("platform_guess")) if config.platform_guess == "unknown" and audit.get("platform_guess") else config.platform_guess),
        coordinate_unit=(str(audit.get("coordinate_unit")) if config.coordinate_unit in {"", "unknown"} and audit.get("coordinate_unit") else config.coordinate_unit),
    )
    source_unit = str(config.coordinate_unit or "native").lower()
    scale = float(config.coordinate_scale) if config.coordinate_scale is not None else None
    has_scale_source = bool(str(config.scale_source).strip())
    physical_calibration = has_scale_source and (source_unit in {"micrometer", "micrometre", "um", "µm"} or (scale is not None and scale > 0))
    if physical_calibration:
        factor = 1.0 if source_unit in {"micrometer", "micrometre", "um", "µm"} else float(scale)
        coords = coords * factor
        config = replace(
            config,
            source_coordinate_unit=source_unit,
            coordinate_unit="micrometer",
            coordinate_scale=factor,
        )
    else:
        config = replace(config, source_coordinate_unit=source_unit)
    if config.radius is not None and not physical_calibration:
        coord_warnings.append("user radius is expressed in native/non-calibrated coordinate units, not calibrated physical distance")
    n_obs = int(getattr(adata, "n_obs", len(coords)))
    if coords.ndim != 2 or coords.shape != (n_obs, 2) or not np.isfinite(coords).all():
        conn, dist = _empty_graph(n_obs)
        qc, qc_warnings, degree_table = graph_qc(
            conn,
            dist,
            None,
            method=requested_method,
            isolated_fraction_warning=config.isolated_fraction_warning,
            largest_component_ratio_warning=config.largest_component_ratio_warning,
            median_degree_warning=config.median_degree_warning,
            long_edge_fraction_warning=config.long_edge_fraction_warning,
            near_complete_density_warning=config.near_complete_density_warning,
        )
        qc["zero_distance_edges_excluded"] = 0
        metadata = _metadata(
            config,
            requested_method,
            requested_method,
            coord_source,
            False,
            qc,
            radius=None,
            lattice_detected=False,
        )
        return GraphBuildResult(requested_method, conn, dist, metadata, qc, coord_warnings + qc_warnings, degree_table)
    effective_method = requested_method
    if requested_method == "radius":
        conn, dist, radius_value, warnings, zero_distance_edges_excluded = _radius_graph(coords, config.radius, config)
        lattice_detected = False
    elif requested_method == "knn":
        conn, dist, warnings, zero_distance_edges_excluded = _knn_graph(coords, config)
        radius_value = config.radius
        lattice_detected = False
    else:
        conn, dist, lattice_detected, effective_method, radius_value, warnings, zero_distance_edges_excluded = _lattice_graph(adata, coords, config)
    spacing = estimate_radius(coords, multiplier=1.0)
    qc, qc_warnings, degree_table = graph_qc(
        conn,
        dist,
        coords,
        method=effective_method,
        typical_spacing=spacing,
        long_edge_factor=config.long_edge_factor,
        isolated_fraction_warning=config.isolated_fraction_warning,
        largest_component_ratio_warning=config.largest_component_ratio_warning,
        median_degree_warning=config.median_degree_warning,
        long_edge_fraction_warning=config.long_edge_fraction_warning,
        near_complete_density_warning=config.near_complete_density_warning,
    )
    qc["zero_distance_edges_excluded"] = int(zero_distance_edges_excluded)
    if zero_distance_edges_excluded:
        warnings.append(
            f"inverse-distance weighting excluded zero-distance edges: {zero_distance_edges_excluded}"
        )
    metadata = _metadata(
        config,
        requested_method,
        effective_method,
        coord_source,
        True,
        qc,
        radius=radius_value,
        lattice_detected=lattice_detected,
    )
    return GraphBuildResult(effective_method, conn, dist, metadata, qc, coord_warnings + warnings + qc_warnings, degree_table)


def _metadata(
    config: GraphBuildConfig,
    requested_method: str,
    effective_method: str,
    coord_source: str,
    coord_valid: bool,
    qc: dict[str, Any],
    *,
    radius: float | None,
    lattice_detected: bool,
) -> dict[str, Any]:
    return json_safe({
        "schema_version": GRAPH_SCHEMA_VERSION,
        "active_graph": effective_method,
        "method": effective_method,
        "requested_method": requested_method,
        "effective_method": effective_method,
        "fallback_used": requested_method != effective_method,
        "radius": radius,
        "graph_method": effective_method,
        "requested_radius": config.radius,
        "effective_radius": radius,
        "radius_unit": config.coordinate_unit or "native",
        "k": int(config.k) if effective_method == "knn" else None,
        "lattice_ring": int(config.lattice_ring) if requested_method == "lattice" else None,
        "lattice_detection_succeeded": bool(lattice_detected),
        "weighting": config.weighting,
        "gaussian_sigma": config.gaussian_sigma,
        "symmetrization": config.symmetrization if effective_method == "knn" else None,
        "coordinate_source": coord_source,
        "source_coordinate_unit": config.source_coordinate_unit or config.coordinate_unit,
        "coordinate_unit": config.coordinate_unit,
        "coordinate_scale": config.coordinate_scale,
        "scale_source": config.scale_source,
        "physical_calibration_available": bool(
            str(config.coordinate_unit).lower() == "micrometer"
            and config.coordinate_scale is not None
            and float(config.coordinate_scale) > 0
            and bool(str(config.scale_source).strip())
        ),
        "platform_guess": config.platform_guess,
        "coordinate_valid": bool(coord_valid),
        "n_nodes": qc.get("n_nodes", 0),
        "n_edges": qc.get("n_edges", 0),
        "mean_degree": qc.get("mean_degree", 0.0),
        "median_degree": qc.get("median_degree", 0.0),
        "isolated_fraction": qc.get("isolated_fraction", 0.0),
        "n_connected_components": qc.get("n_connected_components", 0),
        "largest_component_ratio": qc.get("largest_component_ratio", 0.0),
        "zero_distance_edges_excluded": qc.get("zero_distance_edges_excluded", 0),
        "qc_thresholds": {
            "isolated_spot_fraction": config.isolated_fraction_warning,
            "largest_component_ratio": config.largest_component_ratio_warning,
            "median_degree": config.median_degree_warning,
            "long_edge_fraction": config.long_edge_fraction_warning,
            "near_complete_density": config.near_complete_density_warning,
        },
        "random_seed": int(config.random_seed),
        "parameters": config,
    })


def store_graph(adata, result: GraphBuildResult, *, active: bool = True) -> None:
    """Store graph matrices and SpatialTX metadata in standard AnnData slots."""
    keys = GRAPH_KEYS.get(result.method, GRAPH_KEYS["radius"])
    if not hasattr(adata, "obsp"):
        raise ValueError("AnnData-like object must provide .obsp for graph storage")
    adata.obsp[keys["connectivities"]] = result.connectivities.tocsr()
    adata.obsp[keys["distances"]] = result.distances.tocsr()
    if not hasattr(adata, "uns"):
        raise ValueError("AnnData-like object must provide .uns for graph metadata")
    existing = dict(adata.uns.get("spatialtx_graph", {}))
    graphs = dict(existing.get("graphs", {}))
    graphs[result.method] = result.metadata
    existing.update({
        "schema_version": GRAPH_SCHEMA_VERSION,
        "active_graph": result.method if active else existing.get("active_graph", result.method),
        "graphs": graphs,
        "last_qc": result.qc,
        "warnings": result.warnings,
        "obsp_keys": GRAPH_KEYS,
    })
    adata.uns["spatialtx_graph"] = json_safe(existing)
