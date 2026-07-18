from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


def _configure_plot_font(matplotlib) -> None:
    import platform

    if platform.system() == "Windows":
        matplotlib.rcParams["font.sans-serif"] = ["Malgun Gothic", "DejaVu Sans"]
        matplotlib.rcParams["axes.unicode_minus"] = False


def _save(fig, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    return path


def _plot_edge_indices(connectivities: sparse.spmatrix, max_plot_edges: int, plot_edge_seed: int) -> tuple[np.ndarray, np.ndarray, int]:
    upper = sparse.triu(connectivities, k=1).tocoo()
    total = int(upper.nnz)
    limit = max(1, int(max_plot_edges))
    if total <= limit:
        selected = np.arange(total, dtype=int)
    else:
        selected = np.sort(np.random.default_rng(int(plot_edge_seed)).choice(total, size=limit, replace=False))
    return np.asarray(upper.row[selected], dtype=int), np.asarray(upper.col[selected], dtype=int), total


def plot_graph_qc(
    coords: np.ndarray,
    connectivities: sparse.spmatrix,
    isolated: np.ndarray,
    path: Path,
    title: str,
    *,
    max_plot_edges: int = 50000,
    plot_edge_seed: int = 42,
) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    _configure_plot_font(matplotlib)
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    coords = np.asarray(coords, dtype=float)
    fig, ax = plt.subplots(figsize=(7.5, 7.0))
    edge_i, edge_j, total_edges = _plot_edge_indices(connectivities, max_plot_edges, plot_edge_seed)
    displayed_edges = int(len(edge_i))
    if displayed_edges:
        segments = np.stack([coords[edge_i, :2], coords[edge_j, :2]], axis=1)
        ax.add_collection(LineCollection(segments, colors="#94a3b8", linewidths=0.35, alpha=0.45, zorder=1))
    ax.scatter(coords[:, 0], coords[:, 1], s=12, c="#1d4ed8", linewidths=0, zorder=2, label="connected spot")
    if isolated.any():
        ax.scatter(coords[isolated, 0], coords[isolated, 1], s=24, c="#dc2626", linewidths=0, zorder=3, label="isolated spot")
    ax.set_title(title)
    ax.set_aspect("equal", adjustable="datalim")
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    ax.legend(frameon=False, loc="best")
    fig.text(
        0.02,
        0.01,
        f"Displayed {displayed_edges:,} of {total_edges:,} graph edges; sampling affects visualization only.",
        fontsize=8,
    )
    result = _save(fig, path)
    plt.close(fig)
    metadata_path = Path(path).with_suffix(Path(path).suffix + ".metadata.json")
    metadata_path.write_text(
        json.dumps(
            {
                "total_graph_edges": total_edges,
                "displayed_plot_edges": displayed_edges,
                "deterministic_downsampling_applied": displayed_edges < total_edges,
                "max_plot_edges": int(max_plot_edges),
                "plot_edge_seed": int(plot_edge_seed),
                "statistics_use_full_graph": True,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    return result


def plot_enrichment_heatmap(table: pd.DataFrame, path: Path, title: str) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    _configure_plot_font(matplotlib)
    import matplotlib.pyplot as plt

    if table.empty:
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.text(0.5, 0.5, "No neighborhood enrichment rows", ha="center", va="center")
        ax.axis("off")
        result = _save(fig, path)
        plt.close(fig)
        return result
    labels = list(dict.fromkeys(list(table["label_a"].astype(str)) + list(table["label_b"].astype(str))))
    matrix = np.full((len(labels), len(labels)), np.nan)
    q_matrix = np.full_like(matrix, np.nan, dtype=float)
    index = {label: i for i, label in enumerate(labels)}
    for _, row in table.iterrows():
        a, b = index[str(row["label_a"])], index[str(row["label_b"])]
        matrix[a, b] = matrix[b, a] = float(row.get("enrichment_z", np.nan))
        q_matrix[a, b] = q_matrix[b, a] = float(row.get("fdr_bh", np.nan))
    finite = matrix[np.isfinite(matrix)]
    limit = max(2.0, float(np.nanmax(np.abs(finite))) if len(finite) else 2.0)
    fig, ax = plt.subplots(figsize=(max(6, len(labels) * 0.7), max(5, len(labels) * 0.62)))
    image = ax.imshow(matrix, cmap="coolwarm", vmin=-limit, vmax=limit)
    ax.set_xticks(range(len(labels)), labels=labels, rotation=45, ha="right")
    ax.set_yticks(range(len(labels)), labels=labels)
    ax.set_title(title)
    for y in range(len(labels)):
        for x in range(len(labels)):
            if np.isfinite(q_matrix[y, x]) and q_matrix[y, x] < 0.05:
                ax.text(x, y, "*", ha="center", va="center", color="black", fontsize=12)
    fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04, label="Enrichment z-score")
    fig.text(0.02, 0.01, "Asterisk marks FDR < 0.05; exploratory spatial association, not causal interaction.", fontsize=8)
    fig.tight_layout()
    result = _save(fig, path)
    plt.close(fig)
    return result


def plot_context_map(coords: np.ndarray, values: np.ndarray, path: Path, title: str, cmap: str = "viridis") -> Path:
    import matplotlib
    matplotlib.use("Agg")
    _configure_plot_font(matplotlib)
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.0, 6.7))
    image = ax.scatter(coords[:, 0], coords[:, 1], c=values, s=13, cmap=cmap, linewidths=0)
    ax.set_title(title)
    ax.set_aspect("equal", adjustable="datalim")
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    result = _save(fig, path)
    plt.close(fig)
    return result


def plot_joint_hv_map(coords: np.ndarray, h_high: np.ndarray, v_high: np.ndarray, path: Path, title: str) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    _configure_plot_font(matplotlib)
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    state = np.full(len(coords), "neither", dtype=object)
    state[h_high & ~v_high] = "H-high"
    state[~h_high & v_high] = "V-high"
    state[h_high & v_high] = "H-high + V-high"
    palette = {"neither": "#cbd5e1", "H-high": "#7c3aed", "V-high": "#059669", "H-high + V-high": "#dc2626"}
    colors = [palette[x] for x in state]
    fig, ax = plt.subplots(figsize=(7.0, 6.7))
    ax.scatter(coords[:, 0], coords[:, 1], c=colors, s=13, linewidths=0)
    ax.set_title(title)
    ax.set_aspect("equal", adjustable="datalim")
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    ax.legend(handles=[Patch(color=color, label=label) for label, color in palette.items()], frameon=False, loc="best")
    result = _save(fig, path)
    plt.close(fig)
    return result
