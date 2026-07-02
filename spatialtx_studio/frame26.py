from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.ndimage import binary_dilation, distance_transform_edt, gaussian_filter
from sklearn.neighbors import NearestNeighbors

from .gene_program import select_gene_programs

C_GENES_FIXED = ["CD8A", "CD8B", "NKG7", "PRF1", "GZMB", "IFNG"]
B_GENES_FIXED = ["COL1A1", "COL1A2", "COL3A1", "FN1", "LUM", "DCN"]


@dataclass
class Frame26Config:
    sigma: float = 0.8
    qscale_low: float = 0.05
    qscale_high: float = 0.95
    c_q: float = 0.80
    b_q: float = 0.80
    g_q: float = 0.60
    transition_radius_spots: int = 2
    diffuse_exclude_interface_buffer: int = 1
    eps: float = 1e-9


def frame26_config_from_dict(config: dict | None) -> Frame26Config:
    if not config:
        return Frame26Config()
    frame_cfg = config.get("frame26", {})
    thresholds = config.get("thresholds", {})
    smoothing = config.get("smoothing", {})
    return Frame26Config(
        sigma=float(frame_cfg.get("sigma", smoothing.get("sigma", Frame26Config.sigma))),
        qscale_low=float(frame_cfg.get("qscale_low", thresholds.get("qscale_low", Frame26Config.qscale_low))),
        qscale_high=float(frame_cfg.get("qscale_high", thresholds.get("qscale_high", Frame26Config.qscale_high))),
        c_q=float(frame_cfg.get("c_q", thresholds.get("c_quantile", Frame26Config.c_q))),
        b_q=float(frame_cfg.get("b_q", thresholds.get("b_quantile", Frame26Config.b_q))),
        g_q=float(frame_cfg.get("g_q", thresholds.get("g_quantile", Frame26Config.g_q))),
        transition_radius_spots=int(
            frame_cfg.get(
                "transition_radius_spots",
                config.get("transition_zone", {}).get("radius_spots", Frame26Config.transition_radius_spots),
            )
        ),
        diffuse_exclude_interface_buffer=int(
            frame_cfg.get("diffuse_exclude_interface_buffer", Frame26Config.diffuse_exclude_interface_buffer)
        ),
        eps=float(frame_cfg.get("eps", Frame26Config.eps)),
    )


def _to_numpy(x):
    if hasattr(x, "toarray"):
        x = x.toarray()
    return np.asarray(x)


def get_expr_vector(adata, gene: str):
    if gene not in adata.var_names:
        return None
    return _to_numpy(adata[:, gene].X).reshape(-1)


def mean_gene_score(adata, genes: list[str]) -> np.ndarray:
    values = []
    for gene in genes:
        vector = get_expr_vector(adata, gene)
        if vector is not None:
            values.append(vector)
    if not values:
        return np.zeros(adata.n_obs, dtype=float)
    return np.mean(np.vstack(values), axis=0)


def qscale(values, q_low=0.05, q_high=0.95, eps=1e-9):
    values = np.asarray(values, dtype=float)
    lo = np.nanquantile(values, q_low)
    hi = np.nanquantile(values, q_high)
    if hi - lo < eps:
        return np.zeros_like(values, dtype=float)
    return np.clip((values - lo) / (hi - lo), 0.0, 1.0)


def build_grid(obs, values, xcol="array_col", ycol="array_row"):
    obs = obs.copy()
    xs = obs[xcol].astype(int).to_numpy()
    ys = obs[ycol].astype(int).to_numpy()
    ux = np.sort(np.unique(xs))
    uy = np.sort(np.unique(ys))
    x_to_i = {x: i for i, x in enumerate(ux)}
    y_to_i = {y: i for i, y in enumerate(uy)}
    grid = np.full((len(uy), len(ux)), np.nan, dtype=float)
    mask = np.zeros((len(uy), len(ux)), dtype=bool)
    for k in range(len(obs)):
        i = y_to_i[ys[k]]
        j = x_to_i[xs[k]]
        grid[i, j] = values[k]
        mask[i, j] = True
    return grid, mask, xs, ys, x_to_i, y_to_i


def smooth_on_grid(grid, mask, sigma=0.8):
    filled = np.where(mask, grid, 0.0)
    weights = mask.astype(float)
    numerator = gaussian_filter(filled, sigma=sigma)
    denominator = gaussian_filter(weights, sigma=sigma)
    out = np.full_like(filled, np.nan, dtype=float)
    ok = denominator > 1e-9
    out[ok] = numerator[ok] / denominator[ok]
    return out


def gradient_magnitude(grid):
    filled = np.nan_to_num(grid, nan=0.0)
    gy, gx = np.gradient(filled)
    return np.sqrt(gx**2 + gy**2)


def flatten_grid_to_spots(obs, grid, x_to_i, y_to_i, xcol="array_col", ycol="array_row"):
    xs = obs[xcol].astype(int).to_numpy()
    ys = obs[ycol].astype(int).to_numpy()
    out = np.zeros(len(obs), dtype=float)
    for k in range(len(obs)):
        out[k] = grid[y_to_i[ys[k]], x_to_i[xs[k]]]
    return out


def compute_sign_map(r_values, eps=1e-9):
    r_values = np.asarray(r_values, dtype=float)
    sign = np.zeros_like(r_values, dtype=int)
    sign[r_values > eps] = 1
    sign[r_values < -eps] = -1
    return sign


def quantile_threshold(values, q):
    values = np.asarray(values, dtype=float)
    finite = np.isfinite(values)
    if finite.sum() == 0:
        return np.nan
    return float(np.nanquantile(values[finite], q))


def detect_interface_like_spots(c_values, b_values, g_values, c_q=0.8, b_q=0.8, g_q=0.6):
    c_th = quantile_threshold(c_values, c_q)
    b_th = quantile_threshold(b_values, b_q)
    g_th = quantile_threshold(g_values, g_q)
    mask = (c_values > c_th) & (b_values > b_th) & (g_values > g_th)
    return mask, c_th, b_th, g_th


def spot_mask_to_grid(obs, mask, x_to_i, y_to_i, xcol="array_col", ycol="array_row"):
    xs = obs[xcol].astype(int).to_numpy()
    ys = obs[ycol].astype(int).to_numpy()
    height = max(y_to_i.values()) + 1
    width = max(x_to_i.values()) + 1
    grid = np.zeros((height, width), dtype=bool)
    for k in range(len(obs)):
        grid[y_to_i[ys[k]], x_to_i[xs[k]]] = bool(mask[k])
    return grid


def grid_mask_to_spots(obs, grid, x_to_i, y_to_i, xcol="array_col", ycol="array_row"):
    xs = obs[xcol].astype(int).to_numpy()
    ys = obs[ycol].astype(int).to_numpy()
    out = np.zeros(len(obs), dtype=bool)
    for k in range(len(obs)):
        out[k] = bool(grid[y_to_i[ys[k]], x_to_i[xs[k]]])
    return out


def compute_distance_to_interface(interface_grid):
    if np.sum(interface_grid) == 0:
        return np.full(interface_grid.shape, np.inf, dtype=float)
    return distance_transform_edt(~interface_grid)


def make_transition_shell(distance_grid, radius_spots=2):
    return (distance_grid > 0) & (distance_grid <= radius_spots)


def safe_quantile(values, q):
    values = np.asarray(values, dtype=float)
    finite = np.isfinite(values)
    if finite.sum() == 0:
        return np.nan
    return float(np.nanquantile(values[finite], q))


def robust_scale01(x, denom_floor=1e-9):
    x = float(x)
    if not np.isfinite(x):
        return 0.0
    return x / (x + 1.0 + denom_floor)


def safe_mean(values):
    values = np.asarray(values, dtype=float)
    finite = np.isfinite(values)
    if finite.sum() == 0:
        return 0.0
    return float(np.nanmean(values[finite]))


def compute_sign_fractions(sign_vec):
    sign_vec = np.asarray(sign_vec, dtype=int)
    n = len(sign_vec)
    if n == 0:
        return 0.0, 0.0, 0.0
    pos_fraction = float(np.mean(sign_vec > 0))
    neg_fraction = float(np.mean(sign_vec < 0))
    nonzero = pos_fraction + neg_fraction
    sign_mixing_fraction = 0.0 if nonzero <= 0 else float((2.0 * min(pos_fraction, neg_fraction)) / nonzero)
    return pos_fraction, neg_fraction, sign_mixing_fraction


def build_neighbor_graph(obs, xcol="array_col", ycol="array_row", k=6, radius_scale=1.3):
    coords = obs[[xcol, ycol]].astype(float).to_numpy()
    n = len(coords)
    if n < 2:
        return [], np.nan
    nn = NearestNeighbors(n_neighbors=min(k + 1, n), metric="euclidean")
    nn.fit(coords)
    distances, indices = nn.kneighbors(coords)
    first_nn = distances[:, 1] if distances.shape[1] > 1 else np.full(n, np.nan)
    radius = float(np.nanmedian(first_nn) * radius_scale)
    edges = set()
    for i in range(n):
        for distance, j in zip(distances[i, 1:], indices[i, 1:]):
            if distance <= radius:
                a, b = sorted((int(i), int(j)))
                edges.add((a, b))
    return sorted(edges), radius


def compute_local_adjacency_metrics(obs, sign_vec, xcol="array_col", ycol="array_row", k=6, radius_scale=1.3):
    sign_vec = np.asarray(sign_vec, dtype=int)
    edges, radius = build_neighbor_graph(obs, xcol=xcol, ycol=ycol, k=k, radius_scale=radius_scale)
    adj_total = 0
    adj_cross = 0
    adj_same = 0
    adj_zero = 0
    for i, j in edges:
        a = sign_vec[i]
        b = sign_vec[j]
        adj_total += 1
        if a == 0 or b == 0:
            adj_zero += 1
        elif a * b < 0:
            adj_cross += 1
        else:
            adj_same += 1
    if adj_total == 0:
        return {
            "adj_total_count": 0,
            "adj_cross_count": 0,
            "adj_same_count": 0,
            "adj_zero_count": 0,
            "adj_cross_fraction": 0.0,
            "adj_same_fraction": 0.0,
            "adj_zero_fraction": 0.0,
            "adj_neighbor_radius": float(radius) if np.isfinite(radius) else np.nan,
        }
    return {
        "adj_total_count": int(adj_total),
        "adj_cross_count": int(adj_cross),
        "adj_same_count": int(adj_same),
        "adj_zero_count": int(adj_zero),
        "adj_cross_fraction": float(adj_cross / adj_total),
        "adj_same_fraction": float(adj_same / adj_total),
        "adj_zero_fraction": float(adj_zero / adj_total),
        "adj_neighbor_radius": float(radius),
    }


def compute_fragmentation_metrics(obs, mask_vec, xcol="array_col", ycol="array_row", k=6, radius_scale=1.3, small_k=10):
    mask_vec = np.asarray(mask_vec, dtype=bool)
    idx = np.where(mask_vec)[0]
    total = len(idx)
    if total == 0:
        return {
            "n_diffuse_components": 0,
            "largest_diffuse_component_size": 0,
            "largest_diffuse_component_ratio": 0.0,
            "diffuse_component_entropy": 0.0,
            "small_component_fraction": 0.0,
            "frag_neighbor_radius": np.nan,
        }
    sub_obs = obs.iloc[idx].copy().reset_index(drop=True)
    edges, radius = build_neighbor_graph(sub_obs, xcol=xcol, ycol=ycol, k=k, radius_scale=radius_scale)
    parent = list(range(total))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for i, j in edges:
        union(i, j)
    comp_sizes = {}
    for i in range(total):
        root = find(i)
        comp_sizes[root] = comp_sizes.get(root, 0) + 1
    sizes = np.asarray(list(comp_sizes.values()), dtype=float)
    largest = int(np.max(sizes))
    p = sizes / sizes.sum()
    return {
        "n_diffuse_components": int(len(sizes)),
        "largest_diffuse_component_size": largest,
        "largest_diffuse_component_ratio": float(largest / total),
        "diffuse_component_entropy": float(-(p * np.log(p + 1e-12)).sum()),
        "small_component_fraction": float(sizes[sizes < small_k].sum() / total),
        "frag_neighbor_radius": float(radius),
    }


def assign_regime(interface_spot_count, diffuse_fraction, transition_burden_score):
    if interface_spot_count > 0:
        return "Type_A"
    if (diffuse_fraction >= 0.30) and (transition_burden_score >= 0.16):
        return "Type_B_candidate"
    return "Type_C_candidate"


def compute_regime_confidence(interface_fraction, diffuse_fraction, transition_burden_score, interface_spot_count):
    if interface_spot_count > 0:
        score = 0.6 * interface_fraction + 0.4 * transition_burden_score
    else:
        score = 0.6 * diffuse_fraction + 0.4 * transition_burden_score
    return float(score)


def summarize_one_sample(**kwargs):
    cfg = kwargs["cfg"]
    row = {
        "sample": kwargs["sample"],
        "mode": kwargs["mode"],
        "n_spots": int(kwargs["n_spots"]),
        "sigma": cfg.sigma,
        "qscale_low": cfg.qscale_low,
        "qscale_high": cfg.qscale_high,
        "c_q": cfg.c_q,
        "b_q": cfg.b_q,
        "g_q": cfg.g_q,
        "transition_radius_spots": cfg.transition_radius_spots,
        "c_threshold": float(kwargs["c_th"]) if np.isfinite(kwargs["c_th"]) else np.nan,
        "b_threshold": float(kwargs["b_th"]) if np.isfinite(kwargs["b_th"]) else np.nan,
        "g_threshold": float(kwargs["g_th"]) if np.isfinite(kwargs["g_th"]) else np.nan,
        "interface_spot_count": int(np.sum(kwargs["interface_mask"])),
        "transition_spot_count": int(np.sum(kwargs["transition_mask"])),
        "diffuse_transition_spot_count": int(np.sum(kwargs["diffuse_mask"])),
        "interface_fraction": float(np.mean(kwargs["interface_mask"])),
        "transition_fraction": float(np.mean(kwargs["transition_mask"])),
        "diffuse_transition_fraction": float(np.mean(kwargs["diffuse_mask"])),
        "g_high_fraction": float(kwargs["g_high_fraction"]),
        "sign_mixing_fraction": float(kwargs["sign_mixing_fraction"]),
        "r_variance": float(kwargs["r_variance"]),
        "gr_mean": float(kwargs["gr_mean"]),
        "gr_q90": float(kwargs["gr_q90"]),
        "diffuse_fraction": float(kwargs["diffuse_fraction"]),
        "pos_fraction": float(kwargs["pos_fraction"]),
        "neg_fraction": float(kwargs["neg_fraction"]),
        "boundary_likeness_score": float(kwargs["boundary_likeness_score"]),
        "transition_burden_score": float(kwargs["transition_burden_score"]),
        "regime_label": str(kwargs["regime_label"]),
        "regime_confidence": float(kwargs["regime_confidence"]),
        "warnings": ";".join(kwargs["warnings"]) if kwargs["warnings"] else "",
    }
    row.update(kwargs["adj_metrics"])
    row.update(kwargs["frag_metrics"])
    return pd.DataFrame([row])


def run_frame26(
    input_path: str | Path,
    sample_name: str,
    output_dir: str | Path,
    mode: str = "fixed",
    cfg: Frame26Config | None = None,
    config: dict | None = None,
):
    cfg = cfg or frame26_config_from_dict(config)
    output_dir = Path(output_dir)
    figures_dir = output_dir / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    warnings = []

    adata = sc.read_h5ad(input_path)
    adata.var_names_make_unique()
    required_cols = ["array_row", "array_col"]
    for col in required_cols:
        if col not in adata.obs.columns:
            raise ValueError(f"Missing required obs column: {col}")
    if "in_tissue" in adata.obs.columns:
        adata = adata[adata.obs["in_tissue"].astype(int) == 1].copy()

    if config:
        selected_c, selected_b, selection_note = select_gene_programs(adata, config, mode)
    else:
        selected_c = [g for g in C_GENES_FIXED if g in adata.var_names]
        selected_b = [g for g in B_GENES_FIXED if g in adata.var_names]
        selection_note = "fixed_gene_program"

    c_raw = mean_gene_score(adata, selected_c)
    b_raw = mean_gene_score(adata, selected_b)
    c_score = qscale(c_raw, q_low=cfg.qscale_low, q_high=cfg.qscale_high, eps=cfg.eps)
    b_score = qscale(b_raw, q_low=cfg.qscale_low, q_high=cfg.qscale_high, eps=cfg.eps)
    r_score = c_score - b_score
    sign_map = compute_sign_map(r_score)

    c_grid, tissue_mask, xs, ys, x_to_i, y_to_i = build_grid(adata.obs, c_score)
    b_grid, _, _, _, _, _ = build_grid(adata.obs, b_score)
    r_grid, _, _, _, _, _ = build_grid(adata.obs, r_score)

    c_smooth = smooth_on_grid(c_grid, tissue_mask, sigma=cfg.sigma)
    b_smooth = smooth_on_grid(b_grid, tissue_mask, sigma=cfg.sigma)
    r_smooth = smooth_on_grid(r_grid, tissue_mask, sigma=cfg.sigma)

    g_c_grid = gradient_magnitude(np.nan_to_num(c_smooth, nan=0.0))
    g_b_grid = gradient_magnitude(np.nan_to_num(b_smooth, nan=0.0))
    g_grid = g_c_grid * g_b_grid
    g_r_grid = gradient_magnitude(np.nan_to_num(r_smooth, nan=0.0))

    c_spot = flatten_grid_to_spots(adata.obs, c_smooth, x_to_i, y_to_i)
    b_spot = flatten_grid_to_spots(adata.obs, b_smooth, x_to_i, y_to_i)
    g_spot = flatten_grid_to_spots(adata.obs, g_grid, x_to_i, y_to_i)
    gr_spot = flatten_grid_to_spots(adata.obs, g_r_grid, x_to_i, y_to_i)

    interface_mask, c_th, b_th, g_th = detect_interface_like_spots(
        c_spot, b_spot, g_spot, c_q=cfg.c_q, b_q=cfg.b_q, g_q=cfg.g_q
    )
    if np.sum(interface_mask) == 0:
        warnings.append("no_interface_like_spots_detected")

    interface_grid = spot_mask_to_grid(adata.obs, interface_mask, x_to_i, y_to_i)
    distance_grid = compute_distance_to_interface(interface_grid)
    if not np.isfinite(distance_grid).any() or np.sum(interface_grid) == 0:
        warnings.append("no_finite_distance_to_interface")

    transition_grid = make_transition_shell(distance_grid, radius_spots=cfg.transition_radius_spots) & tissue_mask
    transition_mask = grid_mask_to_spots(adata.obs, transition_grid, x_to_i, y_to_i)
    if np.sum(transition_mask) == 0:
        warnings.append("empty_transition_shell")

    gr_high_th = safe_quantile(gr_spot, cfg.g_q)
    if np.isfinite(gr_high_th):
        high_grad_mask = gr_spot >= gr_high_th
    else:
        high_grad_mask = np.zeros_like(interface_mask, dtype=bool)

    if np.sum(interface_mask) > 0:
        dilated_interface_grid = binary_dilation(interface_grid, iterations=cfg.diffuse_exclude_interface_buffer)
        dilated_interface_mask = grid_mask_to_spots(adata.obs, dilated_interface_grid, x_to_i, y_to_i)
    else:
        dilated_interface_mask = np.zeros_like(interface_mask, dtype=bool)
    diffuse_mask = high_grad_mask & (~interface_mask) & (~dilated_interface_mask)

    g_high_fraction = float(np.mean(high_grad_mask))
    pos_fraction, neg_fraction, sign_mixing_fraction = compute_sign_fractions(sign_map)
    r_variance = float(np.nanvar(r_score))
    scaled_r_variance = robust_scale01(r_variance)
    gr_mean = safe_mean(gr_spot)
    gr_q90 = safe_quantile(gr_spot, 0.90)
    scaled_gr_mean = robust_scale01(gr_mean)
    scaled_gr_q90 = robust_scale01(gr_q90)
    diffuse_fraction = float(np.mean(diffuse_mask))
    interface_fraction = float(np.mean(interface_mask))

    boundary_likeness_score = float(
        0.55 * interface_fraction
        + 0.20 * g_high_fraction
        + 0.15 * scaled_gr_mean
        + 0.10 * scaled_gr_q90
    )
    transition_burden_score = float(
        0.30 * diffuse_fraction
        + 0.20 * scaled_gr_mean
        + 0.20 * scaled_gr_q90
        + 0.15 * sign_mixing_fraction
        + 0.15 * scaled_r_variance
    )
    adj_metrics = compute_local_adjacency_metrics(adata.obs, sign_map)
    frag_metrics = compute_fragmentation_metrics(adata.obs, diffuse_mask)
    regime_label = assign_regime(int(np.sum(interface_mask)), diffuse_fraction, transition_burden_score)
    regime_confidence = compute_regime_confidence(
        interface_fraction, diffuse_fraction, transition_burden_score, int(np.sum(interface_mask))
    )

    adata.obs["frame26_C"] = c_score
    adata.obs["frame26_B"] = b_score
    adata.obs["frame26_R"] = r_score
    adata.obs["frame26_sign"] = sign_map
    adata.obs["frame26_C_smooth"] = c_spot
    adata.obs["frame26_B_smooth"] = b_spot
    adata.obs["frame26_G"] = g_spot
    adata.obs["frame26_GR"] = gr_spot
    adata.obs["frame26_interface_like"] = interface_mask.astype(int)
    adata.obs["frame26_transition_shell"] = transition_mask.astype(int)
    adata.obs["frame26_diffuse_transition"] = diffuse_mask.astype(int)
    adata.obs["frame26_high_grad"] = high_grad_mask.astype(int)

    summary = summarize_one_sample(
        sample=sample_name,
        mode=mode,
        n_spots=adata.n_obs,
        cfg=cfg,
        c_th=c_th,
        b_th=b_th,
        g_th=g_th,
        interface_mask=interface_mask,
        transition_mask=transition_mask,
        diffuse_mask=diffuse_mask,
        g_high_fraction=g_high_fraction,
        sign_mixing_fraction=sign_mixing_fraction,
        r_variance=r_variance,
        boundary_likeness_score=boundary_likeness_score,
        transition_burden_score=transition_burden_score,
        gr_mean=gr_mean,
        gr_q90=gr_q90,
        diffuse_fraction=diffuse_fraction,
        pos_fraction=pos_fraction,
        neg_fraction=neg_fraction,
        regime_label=regime_label,
        regime_confidence=regime_confidence,
        adj_metrics=adj_metrics,
        frag_metrics=frag_metrics,
        warnings=warnings,
    )
    summary["selection_note"] = selection_note
    summary["selected_C"] = ",".join(selected_c)
    summary["selected_B"] = ",".join(selected_b)

    summary.to_csv(output_dir / "metrics.csv", index=False)
    summary.to_csv(output_dir / "frame26_summary.csv", index=False)
    write_frame26_qc(output_dir / "qc_report.csv", summary, warnings)
    pd.DataFrame(
        [{"program": "C", "gene": g} for g in selected_c]
        + [{"program": "B", "gene": g} for g in selected_b]
    ).to_csv(output_dir / "selected_genes.csv", index=False)
    save_frame26_figures(figures_dir, c_grid, b_grid, r_grid, tissue_mask, interface_grid, transition_grid, diffuse_mask, adata.obs, x_to_i, y_to_i)
    with (output_dir / "run_log.txt").open("w", encoding="utf-8") as f:
        for col in summary.columns:
            f.write(f"{col}={summary.iloc[0][col]}\n")
    return summary, adata


def write_frame26_qc(path: Path, summary: pd.DataFrame, warnings: list[str]) -> None:
    row = summary.iloc[0]
    qc = [
        {"check": "analysis", "status": "ok", "detail": "Default C/S balance-field workflow"},
        {"check": "gene_selection", "status": "ok", "detail": str(row["selection_note"])},
        {"check": "regime_label", "status": "ok", "detail": str(row["regime_label"])},
        {"check": "interface_detection", "status": "ok" if row["interface_spot_count"] > 0 else "warning", "detail": f"interface_spot_count={row['interface_spot_count']}"},
        {"check": "diffuse_transition", "status": "ok", "detail": f"diffuse_fraction={row['diffuse_fraction']}"},
        {"check": "warnings", "status": "warning" if warnings else "ok", "detail": ";".join(warnings)},
        {"check": "research_use_only", "status": "notice", "detail": "Not intended for clinical diagnosis or treatment decision-making."},
    ]
    pd.DataFrame(qc).to_csv(path, index=False)


def save_frame26_figures(figures_dir: Path, c_grid, b_grid, r_grid, tissue_mask, interface_grid, transition_grid, diffuse_mask, obs, x_to_i, y_to_i):
    background = np.where(tissue_mask, r_grid, np.nan)
    fig, ax = plt.subplots(figsize=(6, 5), dpi=160)
    image = ax.imshow(background, cmap="coolwarm", interpolation="nearest")
    iy, ix = np.where(interface_grid)
    if len(ix) > 0:
        ax.scatter(ix, iy, s=12, c="black", marker="s", label="Interface-like")
    else:
        ax.text(0.5, 0.5, "No interface-like spots", transform=ax.transAxes, ha="center", va="center", fontsize=11, bbox={"facecolor": "white", "edgecolor": "0.4", "alpha": 0.85})
    ax.set_title("SpatialTX interface-like candidate map")
    ax.set_axis_off()
    fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04, label="R = C - S")
    if len(ix) > 0:
        ax.legend(loc="lower right", fontsize=7)
    fig.tight_layout()
    fig.savefig(figures_dir / "interface_map.png", bbox_inches="tight")
    plt.close(fig)

    diffuse_grid = spot_mask_to_grid(obs, diffuse_mask, x_to_i, y_to_i)
    fig, ax = plt.subplots(figsize=(6, 5), dpi=160)
    ax.imshow(background, cmap="coolwarm", interpolation="nearest", alpha=0.75)
    ty, tx = np.where(transition_grid)
    dy, dx = np.where(diffuse_grid)
    if len(tx) > 0:
        ax.scatter(tx, ty, s=9, c="gold", alpha=0.55, marker="s", label="Transition shell")
    if len(dx) > 0:
        ax.scatter(dx, dy, s=9, c="limegreen", alpha=0.65, marker="s", label="Diffuse transition")
    iy, ix = np.where(interface_grid)
    if len(ix) > 0:
        ax.scatter(ix, iy, s=12, c="black", marker="s", label="Interface-like")
    ax.set_title("SpatialTX transition / diffuse map")
    ax.set_axis_off()
    if len(tx) > 0 or len(dx) > 0 or len(ix) > 0:
        ax.legend(loc="lower right", fontsize=7)
    fig.tight_layout()
    fig.savefig(figures_dir / "transition_zone_map.png", bbox_inches="tight")
    plt.close(fig)
