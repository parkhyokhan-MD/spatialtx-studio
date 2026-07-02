from __future__ import annotations

import numpy as np
from scipy.ndimage import gaussian_filter

from .gene_program import mean_program_score

EPS = 1e-9


def get_spatial_coordinates(adata) -> tuple[np.ndarray, np.ndarray, str]:
    obs = adata.obs
    if "array_row" in obs.columns and "array_col" in obs.columns:
        rows = np.asarray(obs["array_row"]).astype(int)
        cols = np.asarray(obs["array_col"]).astype(int)
        return rows, cols, "array_row/array_col"
    if "spatial" in adata.obsm:
        spatial = np.asarray(adata.obsm["spatial"])
        if spatial.ndim != 2 or spatial.shape[1] < 2:
            raise ValueError("obsm['spatial'] exists but has invalid shape")
        rows = np.round(spatial[:, 1]).astype(int)
        cols = np.round(spatial[:, 0]).astype(int)
        return rows, cols, "obsm['spatial']"
    raise ValueError("No spatial coordinates found: need array_row/array_col or obsm['spatial']")


def qscale(x: np.ndarray, qlow: float = 0.05, qhigh: float = 0.95) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    lo = np.nanquantile(x, qlow)
    hi = np.nanquantile(x, qhigh)
    if not np.isfinite(lo) or not np.isfinite(hi) or hi - lo < EPS:
        return np.zeros_like(x, dtype=float)
    return np.clip((x - lo) / (hi - lo + EPS), 0.0, 1.0)


def make_lattice(values: np.ndarray, rows: np.ndarray, cols: np.ndarray):
    values = np.asarray(values, dtype=float)
    rows = np.asarray(rows).astype(int)
    cols = np.asarray(cols).astype(int)
    r0 = rows.min()
    c0 = cols.min()
    rr = rows - r0
    cc = cols - c0
    height = rr.max() + 1
    width = cc.max() + 1

    sums = np.zeros((height, width), dtype=float)
    counts = np.zeros((height, width), dtype=float)
    mask = np.zeros((height, width), dtype=bool)
    for value, row, col in zip(values, rr, cc):
        sums[row, col] += float(value)
        counts[row, col] += 1.0
        mask[row, col] = True
    grid = np.zeros((height, width), dtype=float)
    grid[mask] = sums[mask] / np.maximum(counts[mask], 1.0)
    return grid, mask, {"r0": int(r0), "c0": int(c0), "height": int(height), "width": int(width)}


def unlattice(grid: np.ndarray, rows: np.ndarray, cols: np.ndarray) -> np.ndarray:
    rows = np.asarray(rows).astype(int)
    cols = np.asarray(cols).astype(int)
    rr = rows - rows.min()
    cc = cols - cols.min()
    out = np.full(len(rows), np.nan, dtype=float)
    for i, (row, col) in enumerate(zip(rr, cc)):
        if 0 <= row < grid.shape[0] and 0 <= col < grid.shape[1]:
            out[i] = grid[row, col]
    return out


def smooth_grid(grid: np.ndarray, mask: np.ndarray, sigma: float) -> np.ndarray:
    fill = np.where(mask, grid, 0.0)
    weights = mask.astype(float)
    smooth_values = gaussian_filter(fill, sigma=sigma)
    smooth_weights = gaussian_filter(weights, sigma=sigma)
    out = np.divide(smooth_values, smooth_weights, out=np.zeros_like(smooth_values), where=smooth_weights > 0)
    out[~mask] = np.nan
    return out


def build_fields(adata, selected_c: list[str], selected_b: list[str], config: dict) -> dict:
    sigma = float(config["smoothing"]["sigma"])
    qlow = float(config["thresholds"]["qscale_low"])
    qhigh = float(config["thresholds"]["qscale_high"])

    c_raw = mean_program_score(adata, selected_c)
    b_raw = mean_program_score(adata, selected_b)
    c_score = qscale(c_raw, qlow, qhigh)
    b_score = qscale(b_raw, qlow, qhigh)
    rows, cols, coord_source = get_spatial_coordinates(adata)

    c_grid, tissue_mask, lattice_meta = make_lattice(c_score, rows, cols)
    b_grid, _, _ = make_lattice(b_score, rows, cols)
    c_smooth_grid = smooth_grid(c_grid, tissue_mask, sigma=sigma)
    b_smooth_grid = smooth_grid(b_grid, tissue_mask, sigma=sigma)

    dcy, dcx = np.gradient(np.nan_to_num(c_smooth_grid, nan=0.0))
    dby, dbx = np.gradient(np.nan_to_num(b_smooth_grid, nan=0.0))
    grad_c_grid = np.sqrt(dcy * dcy + dcx * dcx)
    grad_b_grid = np.sqrt(dby * dby + dbx * dbx)
    g_grid = grad_c_grid * grad_b_grid
    g_grid[~tissue_mask] = np.nan

    return {
        "rows": rows,
        "cols": cols,
        "coord_source": coord_source,
        "lattice_meta": lattice_meta,
        "tissue_mask_grid": tissue_mask,
        "C_raw": c_raw,
        "B_raw": b_raw,
        "C": c_score,
        "B": b_score,
        "C_grid": c_grid,
        "B_grid": b_grid,
        "C_smooth_grid": c_smooth_grid,
        "B_smooth_grid": b_smooth_grid,
        "G_grid": g_grid,
        "C_smooth": unlattice(c_smooth_grid, rows, cols),
        "B_smooth": unlattice(b_smooth_grid, rows, cols),
        "G": unlattice(g_grid, rows, cols),
    }
