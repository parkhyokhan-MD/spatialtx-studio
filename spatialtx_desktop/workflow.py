from __future__ import annotations

import datetime as dt
import json
import math
import shutil
import zipfile
from collections import deque
from pathlib import Path
from typing import Callable, Iterable

import numpy as np
import pandas as pd

from . import __version__


DEFAULT_C_GENES = ["CD8A", "CD8B", "NKG7", "PRF1", "GZMB", "IFNG"]
DEFAULT_S_GENES = ["COL1A1", "COL1A2", "COL3A1", "FN1", "LUM", "DCN"]

C_MARKERS = [
    "CD3D", "CD3E", "CD2", "TRAC", "CD8A", "CD8B", "NKG7", "PRF1", "GZMB",
    "GNLY", "IFNG", "CXCL9", "CXCL10", "CCL5", "CXCL13", "LAG3", "PDCD1",
    "HAVCR2", "TIGIT", "MS4A1", "CD79A", "LYZ", "AIF1", "C1QA", "C1QB",
    "FCGR3A", "ITGAM",
]
S_MARKERS = [
    "COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1", "FN1", "LUM", "DCN",
    "SPARC", "POSTN", "ACTA2", "TAGLN", "MMP2", "MMP11", "THY1", "FAP",
    "PDGFRA", "PDGFRB", "VIM", "ITGA11",
]

Progress = Callable[[str], None]


def parse_gene_text(text: str | Iterable[str]) -> list[str]:
    if not isinstance(text, str):
        values = list(text)
    else:
        values = text.replace(";", ",").replace("\n", ",").split(",")
    result: list[str] = []
    seen: set[str] = set()
    for value in values:
        gene = str(value).strip()
        key = gene.upper()
        if gene and key not in seen:
            result.append(gene)
            seen.add(key)
    return result


def scan_h5ad(folder: str | Path) -> list[Path]:
    root = Path(folder).expanduser()
    if not root.is_dir():
        raise ValueError(f"Input folder does not exist: {root}")
    return sorted((p.resolve() for p in root.rglob("*.h5ad") if p.is_file()), key=lambda p: str(p).lower())


def _read_h5ad(path: str | Path):
    import anndata as ad

    file_path = Path(path)
    if not file_path.is_file() or file_path.suffix.lower() != ".h5ad":
        raise ValueError(f"Not an h5ad file: {file_path}")
    return ad.read_h5ad(file_path)


def _is_count_like(X) -> bool:
    if hasattr(X, "data") and not isinstance(X, np.ndarray):
        data = np.asarray(X.data)
        if data.size == 0:
            return False
        return bool(np.nanmin(data) >= 0 and np.nanmax(data) > 20)
    arr = np.asarray(X)
    return bool(arr.size and np.nanmin(arr) >= 0 and np.nanmax(arr) > 20)


def _dense(X) -> np.ndarray:
    if hasattr(X, "toarray"):
        X = X.toarray()
    return np.asarray(X, dtype=float)


def _zscore_columns(matrix: np.ndarray) -> np.ndarray:
    matrix = np.asarray(matrix, dtype=float)
    mean = np.nanmean(matrix, axis=0)
    std = np.nanstd(matrix, axis=0)
    std[~np.isfinite(std) | (std == 0)] = 1.0
    return (matrix - mean) / std


def _gene_indices(adata, requested: list[str]) -> tuple[list[int], list[str], list[str]]:
    lookup: dict[str, tuple[int, str]] = {}
    for i, gene in enumerate(adata.var_names):
        lookup.setdefault(str(gene).upper(), (i, str(gene)))
    indices: list[int] = []
    present: list[str] = []
    missing: list[str] = []
    for gene in requested:
        hit = lookup.get(gene.upper())
        if hit is None:
            missing.append(gene)
        else:
            indices.append(hit[0])
            present.append(hit[1])
    return indices, present, missing


def _extract_coords(adata) -> tuple[np.ndarray, str]:
    def valid(values) -> np.ndarray | None:
        coords = np.asarray(values, dtype=float)
        if coords.ndim == 2 and coords.shape[0] == adata.n_obs and coords.shape[1] >= 2:
            coords = coords[:, :2]
            if np.isfinite(coords).all():
                return coords
        return None

    if "spatial" in adata.obsm:
        coords = valid(adata.obsm["spatial"])
        if coords is not None:
            return coords, "obsm['spatial']"
    for x, y in [
        ("x", "y"), ("X", "Y"), ("array_col", "array_row"),
        ("pxl_col_in_fullres", "pxl_row_in_fullres"), ("spatial_x", "spatial_y"),
    ]:
        if x in adata.obs.columns and y in adata.obs.columns:
            coords = valid(adata.obs[[x, y]].to_numpy(dtype=float))
            if coords is not None:
                return coords, f"obs[{x},{y}]"
    side = int(math.ceil(math.sqrt(adata.n_obs)))
    return np.column_stack([np.arange(adata.n_obs) % side, np.arange(adata.n_obs) // side]), "fallback_grid"


def _knn(coords: np.ndarray, k: int = 6) -> tuple[list[tuple[int, int]], list[list[int]]]:
    from scipy.spatial import cKDTree

    n = len(coords)
    if np.asarray(coords).ndim != 2 or np.asarray(coords).shape[1] < 2 or not np.isfinite(coords).all():
        raise ValueError("Spatial coordinates must be a finite n-by-2 matrix.")
    if n < 2:
        return [], [[] for _ in range(n)]
    _, neighbors = cKDTree(coords).query(coords, k=min(k + 1, n))
    if neighbors.ndim == 1:
        neighbors = neighbors[:, None]
    edges: set[tuple[int, int]] = set()
    for i, row in enumerate(neighbors):
        for j in row:
            j = int(j)
            if i != j:
                edges.add((min(i, j), max(i, j)))
    adj = [[] for _ in range(n)]
    for a, b in sorted(edges):
        adj[a].append(b)
        adj[b].append(a)
    return sorted(edges), adj


def _gradient(values: np.ndarray, adj: list[list[int]]) -> np.ndarray:
    return np.asarray([
        float(np.mean(np.abs(values[i] - values[row]))) if row else 0.0
        for i, row in enumerate(adj)
    ])


def _components(mask: np.ndarray, adj: list[list[int]]) -> list[int]:
    seen = np.zeros(len(mask), dtype=bool)
    sizes: list[int] = []
    for start in np.flatnonzero(mask):
        if seen[start]:
            continue
        queue = deque([int(start)])
        seen[start] = True
        size = 0
        while queue:
            node = queue.popleft()
            size += 1
            for other in adj[node]:
                if mask[other] and not seen[other]:
                    seen[other] = True
                    queue.append(other)
        sizes.append(size)
    return sorted(sizes, reverse=True)


def _adjacency_metrics(R: np.ndarray, edges: list[tuple[int, int]]) -> tuple[float, float, float, float]:
    if not edges:
        return (float("nan"),) * 4
    neutral = float(np.nanquantile(np.abs(R), 0.20))
    sign = np.zeros(len(R), dtype=int)
    sign[R > neutral] = 1
    sign[R < -neutral] = -1
    same = zero = opposite = crossing = 0
    for a, b in edges:
        sa, sb = sign[a], sign[b]
        zero += int(sa == 0 or sb == 0)
        same += int(sa != 0 and sa == sb)
        opposite += int(sa * sb == -1)
        crossing += int((R[a] <= 0 < R[b]) or (R[b] <= 0 < R[a]))
    total = len(edges)
    return same / total, zero / total, opposite / total, crossing / total


def _high_quantile(values: np.ndarray, quantile: float) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    finite = np.isfinite(values)
    result = np.zeros(len(values), dtype=bool)
    if finite.sum() < 2:
        return result
    dynamic_range = float(np.nanmax(values[finite]) - np.nanmin(values[finite]))
    if dynamic_range <= np.finfo(float).eps * 100:
        return result
    threshold = float(np.nanquantile(values[finite], quantile))
    result[finite] = values[finite] >= threshold
    return result


def score_h5ad(
    path: str | Path,
    c_genes: list[str],
    s_genes: list[str],
    c_q: float = 0.80,
    s_q: float = 0.80,
    g_q: float = 0.60,
):
    adata = _read_h5ad(path)
    for name, value in (("C", c_q), ("S", s_q), ("G", g_q)):
        if not 0 < float(value) < 1:
            raise ValueError(f"{name} quantile must be between 0 and 1.")
    if adata.n_obs < 2:
        raise ValueError("Spatial scoring requires at least two spots/observations.")
    c_genes, s_genes = parse_gene_text(c_genes), parse_gene_text(s_genes)
    c_idx, c_present, c_missing = _gene_indices(adata, c_genes)
    s_idx, s_present, s_missing = _gene_indices(adata, s_genes)
    if not c_idx or not s_idx:
        raise ValueError(f"Insufficient C/S genes: C present={len(c_idx)}, S present={len(s_idx)}")
    all_idx = c_idx + s_idx
    expression = _dense(adata.X[:, all_idx])
    if not np.isfinite(expression).any():
        raise ValueError("Selected C/S genes contain no finite expression values.")
    count_like = _is_count_like(adata.X)
    if count_like:
        expression = np.log1p(expression)
    z = _zscore_columns(expression)
    C = np.nanmean(z[:, :len(c_idx)], axis=1)
    S = np.nanmean(z[:, len(c_idx):], axis=1)
    R = C - S
    if not (np.isfinite(C).all() and np.isfinite(S).all() and np.isfinite(R).all()):
        raise ValueError("Selected C/S genes produce non-finite program scores.")
    coords, coord_source = _extract_coords(adata)
    edges, adj = _knn(coords)
    G = _gradient(R, adj)
    high_c = _high_quantile(C, c_q)
    high_s = _high_quantile(S, s_q)
    high_g = _high_quantile(G, g_q)
    interface = high_c & high_s & high_g
    diffuse = high_g & (high_c | high_s) & ~interface

    interface_sizes = _components(interface, adj)
    diffuse_sizes = _components(diffuse, adj)
    interface_spots, diffuse_spots = int(interface.sum()), int(diffuse.sum())
    interface_largest = interface_sizes[0] / interface_spots if interface_spots else 0.0
    diffuse_largest = diffuse_sizes[0] / diffuse_spots if diffuse_spots else 0.0
    small_fraction = (sum(size <= 3 for size in diffuse_sizes) / len(diffuse_sizes)) if diffuse_sizes else 0.0
    interface_fraction = float(np.mean(interface))
    diffuse_fraction = float(np.mean(diffuse))
    adj_same, adj_zero, adj_opposite, r_crossing = _adjacency_metrics(R, edges)
    burden = 0.50 * diffuse_fraction + 0.25 * float(np.nanquantile(G, 0.90)) + 0.25 * (1 - diffuse_largest if diffuse_spots else 0)
    if interface_fraction >= 0.01:
        regime, localized, transition = "Type_A_candidate", "Localized interface-like candidate detected", "Localized interface-like transition pattern"
    elif diffuse_fraction >= 0.05:
        regime, localized, transition = "Type_B_candidate", "Localized interface-like signal not prominent; diffuse transition burden present", "Diffuse transition-zone organization"
    else:
        regime, localized, transition = "Type_C_candidate", "Transition-poor / flat", "No clear transition-zone organization"
    public_pattern = ""
    if regime == "Type_B_candidate":
        if adj_zero >= 0.55 and small_fraction >= 0.20 and diffuse_largest <= 0.35:
            public_pattern = "Fragmented diffuse transition"
        elif adj_same >= 0.80 and adj_zero <= 0.20:
            public_pattern = "Continuous diffuse transition"
        else:
            public_pattern = "Weakly organized diffuse transition"
    c_coverage = len(c_present) / max(1, len(c_genes))
    s_coverage = len(s_present) / max(1, len(s_genes))
    qc_notes: list[str] = []
    if min(c_coverage, s_coverage) < .5:
        qc_notes.append("low_gene_program_coverage")
    elif min(c_coverage, s_coverage) < .8:
        qc_notes.append("partial_gene_program_coverage")
    if coord_source == "fallback_grid":
        qc_notes.append("spatial_coordinates_missing_or_nonfinite")
    case_insensitive_unique = len({str(gene).upper() for gene in adata.var_names}) == adata.n_vars
    if not case_insensitive_unique:
        qc_notes.append("duplicate_feature_names")
    if adata.n_obs < 10:
        qc_notes.append("very_small_spot_count")
    overlap = sorted({gene.upper() for gene in c_present} & {gene.upper() for gene in s_present})
    if overlap:
        qc_notes.append("overlapping_C_S_genes")
    qc_flag = "PASS" if not qc_notes else "FAIL" if "low_gene_program_coverage" in qc_notes else "WARN"
    metrics = {
        "sample": Path(path).stem, "source_h5ad": str(Path(path).resolve()), "status": "ok",
        "n_spots": int(adata.n_obs), "n_genes": int(adata.n_vars), "coordinate_source": coord_source,
        "expression_transform": "log1p_count_like" if count_like else "existing_processed_scale",
        "spatial_neighbors_k": 6,
        "regime_label": regime, "localized_interface_call": localized, "transition_zone_call": transition,
        "public_transition_pattern": public_pattern, "interface_fraction": interface_fraction,
        "interface_spots": interface_spots, "n_interface_components": len(interface_sizes),
        "largest_interface_component_ratio": interface_largest,
        "interface_fragmentation_index": len(interface_sizes) / interface_spots if interface_spots else 0.0,
        "interface_coherence_score": interface_fraction * interface_largest,
        "diffuse_fraction": diffuse_fraction, "diffuse_spots": diffuse_spots,
        "n_diffuse_components": len(diffuse_sizes), "largest_diffuse_component_ratio": diffuse_largest,
        "small_component_fraction": small_fraction, "diffuse_coherence_score": diffuse_fraction * diffuse_largest,
        "transition_burden_score": burden, "adj_same_fraction": adj_same, "adj_zero_fraction": adj_zero,
        "adj_opposite_fraction": adj_opposite, "R_crossing_fraction": r_crossing,
        "C_mean": float(np.nanmean(C)), "S_mean": float(np.nanmean(S)), "R_mean": float(np.nanmean(R)),
        "R_sd": float(np.nanstd(R)), "R_dynamic_range": float(np.nanquantile(R, .9) - np.nanquantile(R, .1)),
        "G_mean": float(np.nanmean(G)), "G_q90": float(np.nanquantile(G, .9)),
        "C_gene_set": ";".join(c_genes), "S_gene_set": ";".join(s_genes),
        "C_genes_present": ";".join(c_present), "S_genes_present": ";".join(s_present),
        "C_genes_missing": ";".join(c_missing), "S_genes_missing": ";".join(s_missing),
        "C_gene_coverage": c_coverage, "S_gene_coverage": s_coverage,
        "QC_flag": qc_flag, "QC_notes": ";".join(qc_notes),
        "C_S_overlap_genes": ";".join(overlap),
        "unique_feature_names": case_insensitive_unique,
        "interface_c_q": c_q, "interface_s_q": s_q, "interface_g_q": g_q,
    }
    fields = {"coords": coords, "C": C, "S": S, "R": R, "G": G, "interface": interface, "diffuse": diffuse}
    return metrics, fields


def save_spatial_map(path: str | Path, metrics: dict, fields: dict) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    coords = fields["coords"]
    fig, axes = plt.subplots(2, 3, figsize=(13, 8))
    panels = [
        (fields["C"], "C(x): immune-side score", "viridis", False),
        (fields["S"], "S(x): stromal-side score", "magma", False),
        (fields["R"], "R(x)=C(x)-S(x)", "coolwarm", False),
        (fields["G"], "G(x): local balance gradient", "plasma", False),
        (fields["interface"], "Localized interface-like candidates", "Greys", True),
        (fields["diffuse"], "Diffuse transition-zone burden", "YlOrRd", True),
    ]
    for ax, (values, title, cmap, binary) in zip(axes.ravel(), panels):
        image = ax.scatter(coords[:, 0], coords[:, 1], c=np.asarray(values, dtype=float), s=10, cmap=cmap,
                           vmin=0 if binary else None, vmax=1 if binary else None, linewidths=0)
        ax.set_title(title, fontsize=10)
        ax.set_aspect("equal", adjustable="datalim")
        ax.invert_yaxis()
        ax.set_xticks([]); ax.set_yticks([])
        fig.colorbar(image, ax=ax, fraction=.046, pad=.03)
    fig.suptitle(
        f"{metrics['sample']} | {metrics['regime_label']}\n"
        f"interface={metrics['interface_fraction']:.3f}, coherence={metrics['interface_coherence_score']:.3f}, "
        f"diffuse={metrics['diffuse_fraction']:.3f}, burden={metrics['transition_burden_score']:.3f}", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, .93])
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return output


def run_batch(paths: list[Path], output_root: str | Path, c_genes: list[str], s_genes: list[str],
              progress: Progress | None = None, c_q: float = .80, s_q: float = .80, g_q: float = .60) -> tuple[Path, pd.DataFrame]:
    if not paths:
        raise ValueError("Select at least one h5ad sample.")
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    run_dir = Path(output_root).expanduser().resolve() / f"spatialtx_run_{stamp}"
    run_dir.mkdir(parents=True, exist_ok=False)
    created = dt.datetime.now().isoformat(timespec="seconds")
    (run_dir / "run_config.json").write_text(json.dumps({
        "application": "SpatialTX Studio Desktop", "version": __version__, "created": created,
        "input_files": [str(Path(path).resolve()) for path in paths],
        "C_gene_program": list(c_genes), "S_gene_program": list(s_genes),
        "interface_quantiles": {"C": c_q, "S": s_q, "G": g_q},
        "spatial_neighbors_k": 6,
    }, indent=2), encoding="utf-8")
    rows: list[dict] = []
    sample_counts: dict[str, int] = {}
    for number, path in enumerate(paths, 1):
        if progress:
            progress(f"[{number}/{len(paths)}] Scoring {path.name}")
        base_name = path.stem
        sample_counts[base_name] = sample_counts.get(base_name, 0) + 1
        sample_name = base_name if sample_counts[base_name] == 1 else f"{base_name}_{sample_counts[base_name]}"
        sample_dir = run_dir / sample_name
        sample_dir.mkdir(parents=True, exist_ok=True)
        try:
            metrics, fields = score_h5ad(path, c_genes, s_genes, c_q, s_q, g_q)
            metrics["source_sample_name"] = base_name
            metrics["sample"] = sample_name
            png = save_spatial_map(sample_dir / f"{sample_name}_spatialtx_maps.png", metrics, fields)
            metrics["spatial_map_png"] = str(png)
            pd.DataFrame([metrics]).to_csv(sample_dir / "metrics.csv", index=False)
            pd.DataFrame(
                [{"program": "C", "gene": g} for g in c_genes] + [{"program": "S", "gene": g} for g in s_genes]
            ).to_csv(sample_dir / "selected_genes.csv", index=False)
            rows.append(metrics)
        except Exception as exc:
            rows.append({"sample": sample_name, "source_sample_name": base_name, "source_h5ad": str(path), "status": f"error: {exc}"})
            if progress:
                progress(f"  Error: {exc}")
    summary = pd.DataFrame(rows)
    summary.to_csv(run_dir / "spatialtx_summary.csv", index=False)
    successful = sum(str(row.get("status", "")) == "ok" for row in rows)
    (run_dir / "RUN_INFO.txt").write_text(
        f"SpatialTX Studio Desktop v{__version__}\nResearch prototype; not for clinical decision-making.\n"
        f"Created: {created}\nSamples: {len(paths)}\nSuccessful: {successful}\nFailed: {len(paths) - successful}\n",
        encoding="utf-8",
    )
    if progress:
        progress(f"Completed: {run_dir}")
    return run_dir, summary


def _minmax(values: np.ndarray) -> np.ndarray:
    values = np.nan_to_num(np.asarray(values, dtype=float))
    low, high = float(np.min(values)), float(np.max(values))
    return (values - low) / (high - low) if high > low else np.zeros_like(values)


def _corr(a: np.ndarray, b: np.ndarray) -> float:
    if np.nanstd(a) == 0 or np.nanstd(b) == 0:
        return 0.0
    value = np.corrcoef(a, b)[0, 1]
    return float(value) if np.isfinite(value) else 0.0


def _energy(x: np.ndarray, diagonal: np.ndarray, pair: np.ndarray) -> float:
    selected = np.flatnonzero(x)
    return float(diagonal @ x + sum(pair[i, j] for n, i in enumerate(selected) for j in selected[n + 1:]))


def _anneal(diagonal: np.ndarray, pair: np.ndarray, k: int, iterations: int, seed: int = 20260624):
    rng = np.random.default_rng(seed)
    k = max(1, min(k, len(diagonal)))
    x = np.zeros(len(diagonal), dtype=int)
    x[np.argsort(diagonal)[:k]] = 1
    current = best = _energy(x, diagonal, pair)
    best_x = x.copy()
    for step in range(iterations):
        chosen, free = np.flatnonzero(x), np.flatnonzero(1 - x)
        if not len(free):
            break
        proposal = x.copy()
        proposal[rng.choice(chosen)] = 0
        proposal[rng.choice(free)] = 1
        candidate = _energy(proposal, diagonal, pair)
        temperature = max(.01, 1 - step / max(1, iterations))
        if candidate < current or rng.random() < np.exp(-(candidate - current) / temperature):
            x, current = proposal, candidate
            if current < best:
                best_x, best = x.copy(), current
    return best_x, best


def _technical(gene: str) -> bool:
    key = gene.upper()
    return key.startswith(("MT-", "RPL", "RPS", "MALAT1"))


def optimize_genes(path: str | Path, side: str, c_genes: list[str], s_genes: list[str], k: int = 8,
                   pool_size: int = 40, iterations: int = 300,
                   candidate_genes: list[str] | None = None, seed: int = 20260624) -> tuple[list[str], pd.DataFrame, dict]:
    side = side.upper()
    if side not in {"C", "S"}:
        raise ValueError("Optimizer side must be C or S.")
    if iterations < 1:
        raise ValueError("Optimizer iterations must be at least 1.")
    baseline, fields = score_h5ad(path, c_genes, s_genes)
    adata = _read_h5ad(path)
    genes = [str(g) for g in adata.var_names]
    lookup: dict[str, str] = {}
    for gene in genes:
        lookup.setdefault(gene.upper(), gene)
    markers = C_MARKERS if side == "C" else S_MARKERS
    current = c_genes if side == "C" else s_genes
    opposite = {g.upper() for g in (s_genes if side == "C" else c_genes)}
    candidates: list[str] = []
    seed_pool = (list(candidate_genes) if candidate_genes else markers) + list(current)
    for gene in seed_pool:
        actual = lookup.get(gene.upper())
        if actual and actual.upper() not in opposite and actual not in candidates and not _technical(actual):
            candidates.append(actual)
    # Add high-variance genes only when marker coverage leaves room in the bounded pool.
    if not candidate_genes and len(candidates) < pool_size:
        X = adata.X
        if hasattr(X, "multiply"):
            mean = np.asarray(X.mean(axis=0)).ravel()
            variance = np.asarray(X.multiply(X).mean(axis=0)).ravel() - mean * mean
        else:
            variance = np.nanvar(np.asarray(X, dtype=float), axis=0)
        for index in np.argsort(-np.nan_to_num(variance)):
            gene = genes[int(index)]
            if gene.upper() not in opposite and gene not in candidates and not _technical(gene):
                candidates.append(gene)
                if len(candidates) >= pool_size:
                    break
    candidates = candidates[:pool_size]
    if len(candidates) < 2:
        raise ValueError("Candidate pool is too small for optimization.")
    if k < 1:
        raise ValueError("Selected gene count (k) must be at least 1.")
    if len(candidates) < k:
        raise ValueError(f"Cannot select exactly k={k} genes from {len(candidates)} available candidates.")
    indices = [genes.index(g) for g in candidates]
    expression = _dense(adata.X[:, indices])
    if _is_count_like(adata.X):
        expression = np.log1p(expression)
    Z = _zscore_columns(expression)
    target = fields[side]
    opposite_target = fields["S" if side == "C" else "C"]
    balance = fields["R"] if side == "C" else -fields["R"]
    target_reward = np.asarray([max(0, _corr(Z[:, i], target)) for i in range(Z.shape[1])])
    balance_reward = np.asarray([max(0, _corr(Z[:, i], balance)) for i in range(Z.shape[1])])
    opposite_penalty = np.asarray([max(0, _corr(Z[:, i], opposite_target)) for i in range(Z.shape[1])])
    gradient_reward = np.asarray([max(0, _corr(Z[:, i], fields["G"])) for i in range(Z.shape[1])])
    interface_mask = np.asarray(fields["interface"], dtype=bool)
    diffuse_mask = np.asarray(fields["diffuse"], dtype=bool)
    if interface_mask.sum() >= 2 and (~interface_mask).sum() >= 2:
        interface_reward = np.nanmean(Z[interface_mask], axis=0) - np.nanmean(Z[~interface_mask], axis=0)
    elif diffuse_mask.sum() >= 2 and (~diffuse_mask).sum() >= 2:
        interface_reward = np.nanmean(Z[diffuse_mask], axis=0) - np.nanmean(Z[~diffuse_mask], axis=0)
    else:
        interface_reward = np.zeros(Z.shape[1], dtype=float)
    r_low, r_high = np.nanquantile(fields["R"], [.2, .8])
    low, high = fields["R"] <= r_low, fields["R"] >= r_high
    enrichment = (np.nanmean(Z[high], axis=0) - np.nanmean(Z[low], axis=0)) if side == "C" else (np.nanmean(Z[low], axis=0) - np.nanmean(Z[high], axis=0))
    detection = np.mean(expression > 0, axis=0)
    reward = (
        .25 * _minmax(target_reward) + .25 * _minmax(balance_reward)
        + .15 * _minmax(np.maximum(enrichment, 0)) + .10 * _minmax(gradient_reward)
        + .10 * _minmax(np.maximum(interface_reward, 0))
        + .075 * _minmax(np.nanvar(Z, axis=0)) + .075 * _minmax(detection)
    )
    with np.errstate(invalid="ignore", divide="ignore"):
        redundancy = np.nan_to_num(np.abs(np.corrcoef(Z.T)))
    np.fill_diagonal(redundancy, 0)
    diagonal = -reward + .25 * (1 - detection) + .20 * _minmax(opposite_penalty)
    selected_mask, energy = _anneal(diagonal, .45 * redundancy, k, iterations, seed=seed)
    selected = [candidates[i] for i in np.flatnonzero(selected_mask)]
    opt_c, opt_s = (selected, s_genes) if side == "C" else (c_genes, selected)
    optimized, _ = score_h5ad(path, opt_c, opt_s)
    detail = pd.DataFrame({
        "gene": candidates, "selected": selected_mask, "side": side, "qubo_reward": reward,
        "target_corr_reward": target_reward, "balance_corr_reward": balance_reward,
        "side_enrichment_reward": enrichment, "gradient_corr_reward": gradient_reward,
        "interface_enrichment_reward": interface_reward,
        "opposite_corr_penalty": opposite_penalty, "detection_fraction": detection, "q_diag": diagonal,
    }).sort_values(["selected", "qubo_reward"], ascending=[False, False])
    summary = {
        "sample": Path(path).stem, "side": side, "selected_genes": ";".join(selected),
        "candidate_pool_size": len(candidates), "requested_k": k, "selected_gene_count": len(selected),
        "iterations": iterations, "random_seed": seed, "qubo_energy": energy,
        "baseline_regime_label": baseline["regime_label"], "optimized_regime_label": optimized["regime_label"],
        "baseline_interface_fraction": baseline["interface_fraction"], "optimized_interface_fraction": optimized["interface_fraction"],
        "delta_interface_fraction": optimized["interface_fraction"] - baseline["interface_fraction"],
        "baseline_R_dynamic_range": baseline["R_dynamic_range"], "optimized_R_dynamic_range": optimized["R_dynamic_range"],
    }
    return selected, detail, summary


def save_optimizer_results(output_root: str | Path, sample: str, side: str, detail: pd.DataFrame, summary: dict) -> Path:
    folder = Path(output_root) / "optimizer"
    folder.mkdir(parents=True, exist_ok=True)
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    detail.to_csv(folder / f"{sample}_{side}_qubo_detail_{stamp}.csv", index=False)
    pd.DataFrame([summary]).to_csv(folder / f"{sample}_{side}_qubo_summary_{stamp}.csv", index=False)
    return folder


def export_result(run_dir: str | Path, destination: str | Path, as_zip: bool) -> Path:
    source = Path(run_dir).resolve()
    destination = Path(destination).expanduser().resolve()
    if not source.is_dir():
        raise ValueError("No completed result folder is available to export.")
    if as_zip:
        target = destination if destination.suffix.lower() == ".zip" else destination.with_suffix(".zip")
        if target == source or source in target.parents:
            raise ValueError("Export ZIP must be saved outside the result folder being archived.")
        target.parent.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(target, "w", zipfile.ZIP_DEFLATED) as archive:
            for file in source.rglob("*"):
                if file.is_file():
                    archive.write(file, Path(source.name) / file.relative_to(source))
        return target
    target = destination / source.name if destination.is_dir() else destination
    if target == source or source in target.parents:
        raise ValueError("Export folder must be outside the source result folder.")
    if target.exists():
        raise FileExistsError(f"Export target already exists: {target}")
    shutil.copytree(source, target)
    return target
