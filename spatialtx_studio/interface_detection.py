from __future__ import annotations

import numpy as np


def detect_interface(fields: dict, config: dict) -> dict:
    c_grid = fields["C_smooth_grid"]
    b_grid = fields["B_smooth_grid"]
    g_grid = fields["G_grid"]
    tissue_mask = fields["tissue_mask_grid"]
    seed_q = float(config["thresholds"]["c_quantile"])
    b_q = float(config["thresholds"]["b_quantile"])
    grad_q = float(config["thresholds"]["g_quantile"])

    valid = tissue_mask & np.isfinite(c_grid) & np.isfinite(b_grid) & np.isfinite(g_grid)
    if valid.sum() == 0:
        raise ValueError("No valid tissue grid points for interface detection")

    c_threshold = float(np.nanquantile(c_grid[valid], seed_q))
    b_threshold = float(np.nanquantile(b_grid[valid], b_q))
    g_threshold = float(np.nanquantile(g_grid[valid], grad_q))
    interface_grid = (c_grid > c_threshold) & (b_grid > b_threshold) & (g_grid > g_threshold) & valid

    interface_spot = np.zeros(len(fields["rows"]), dtype=bool)
    r0 = fields["rows"].min()
    c0 = fields["cols"].min()
    for i, (row, col) in enumerate(zip(fields["rows"] - r0, fields["cols"] - c0)):
        interface_spot[i] = bool(interface_grid[row, col])

    return {
        "c_threshold": c_threshold,
        "b_threshold": b_threshold,
        "g_threshold": g_threshold,
        "interface_grid": interface_grid,
        "interface_spot": interface_spot,
        "interface_detected": bool(interface_grid.sum() > 0),
    }
