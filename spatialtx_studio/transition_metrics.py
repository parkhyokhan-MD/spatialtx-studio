from __future__ import annotations

import numpy as np
from scipy.stats import linregress, pearsonr

EPS = 1e-9


def safe_pearson(x: np.ndarray, y: np.ndarray) -> float:
    valid = np.isfinite(x) & np.isfinite(y)
    if valid.sum() < 3:
        return np.nan
    if np.nanstd(x[valid]) < EPS or np.nanstd(y[valid]) < EPS:
        return np.nan
    return float(pearsonr(x[valid], y[valid])[0])


def estimate_alpha_from_decay(r_field: np.ndarray, distance: np.ndarray, min_points: int = 30):
    r_field = np.asarray(r_field, dtype=float)
    distance = np.asarray(distance, dtype=float)
    valid = np.isfinite(r_field) & np.isfinite(distance) & (distance > 0)
    if valid.sum() < min_points:
        valid = np.isfinite(r_field) & np.isfinite(distance) & (distance >= 0)
    if valid.sum() < min_points:
        return np.nan, np.nan, np.nan, np.nan, np.nan

    x = distance[valid]
    y = np.log(np.abs(r_field[valid]) + EPS)
    if np.nanstd(x) < EPS or np.nanstd(y) < EPS:
        return np.nan, np.nan, np.nan, np.nan, np.nan
    slope, intercept, r_value, _p_value, _std_err = linregress(x, y)
    if not np.isfinite(slope):
        return np.nan, np.nan, np.nan, np.nan, np.nan
    alpha = float(max(-slope, EPS))
    lambda_value = float(1.0 / alpha)
    return alpha, lambda_value, float(slope), float(intercept), float(r_value)


def compute_integrals(r_field: np.ndarray, distance: np.ndarray, alpha: float):
    valid = np.isfinite(r_field) & np.isfinite(distance)
    if valid.sum() == 0 or not np.isfinite(alpha):
        return np.nan, np.nan, np.nan
    r_valid = r_field[valid]
    weights = np.exp(-alpha * distance[valid])
    i_total = float(np.sum(r_valid * weights))
    i_bal = float(np.sum(r_valid[r_valid < 0] * weights[r_valid < 0]))
    i_str = float(np.sum(r_valid[r_valid > 0] * weights[r_valid > 0]))
    return i_total, i_bal, i_str


def compute_delta_interface(c_field: np.ndarray, b_field: np.ndarray, interface_spot: np.ndarray) -> float:
    mask = np.asarray(interface_spot, dtype=bool)
    if mask.sum() == 0 or (~mask).sum() == 0:
        return np.nan
    c_int = np.nanmedian(c_field[mask])
    b_int = np.nanmedian(b_field[mask])
    c_non = np.nanmedian(c_field[~mask])
    b_non = np.nanmedian(b_field[~mask])
    return float((b_int - c_int) - (b_non - c_non))


def compute_metrics(sample: str, gene_mode: str, adata, selected_c: list[str], selected_b: list[str], fields: dict, interface: dict, transition: dict) -> dict:
    c_smooth = fields["C_smooth"]
    b_smooth = fields["B_smooth"]
    r_field = c_smooth - b_smooth
    distance = transition["distance_spot"]
    alpha, lambda_value, slope, intercept, fit_r = estimate_alpha_from_decay(r_field, distance)
    i_total, i_bal, i_str = compute_integrals(r_field, distance, alpha)
    interface_fraction = float(np.mean(interface["interface_spot"])) if len(interface["interface_spot"]) else 0.0

    return {
        "sample": sample,
        "gene_mode": gene_mode,
        "n_spots": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "selected_C": ",".join(selected_c),
        "selected_B": ",".join(selected_b),
        "interface_detected": bool(interface["interface_detected"]),
        "interface_fraction": interface_fraction,
        "alpha": alpha,
        "lambda": lambda_value,
        "I_bal": i_bal,
        "I_str": i_str,
        "delta_interface": compute_delta_interface(c_smooth, b_smooth, interface["interface_spot"]),
        "C_B_correlation": safe_pearson(c_smooth, b_smooth),
        "qc_status": "ok" if interface["interface_detected"] else "no_interface",
        "I_total": i_total,
        "alpha_fit_slope": slope,
        "alpha_fit_intercept": intercept,
        "alpha_fit_r": fit_r,
        "c_threshold": interface["c_threshold"],
        "b_threshold": interface["b_threshold"],
        "g_threshold": interface["g_threshold"],
    }
