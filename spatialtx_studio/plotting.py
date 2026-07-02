from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def save_interface_map(path: str | Path, fields: dict, interface: dict) -> None:
    fig, ax = plt.subplots(figsize=(6, 5), dpi=160)
    c_grid = fields["C_smooth_grid"]
    b_grid = fields["B_smooth_grid"]
    background = c_grid - b_grid
    image = ax.imshow(background, cmap="coolwarm", interpolation="nearest")
    yy, xx = np.where(interface["interface_grid"])
    if len(xx) > 0:
        ax.scatter(xx, yy, s=10, c="black", marker="s", label="Interface")
    else:
        ax.text(
            0.5,
            0.5,
            "No interface detected",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=11,
            bbox={"facecolor": "white", "edgecolor": "0.4", "alpha": 0.85},
        )
    ax.set_title("SpatialTX interface-like candidate map")
    ax.set_axis_off()
    fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04, label="R = C - S")
    if len(xx) > 0:
        ax.legend(loc="lower right", fontsize=7)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def save_transition_zone_map(path: str | Path, fields: dict, interface: dict, transition: dict) -> None:
    fig, ax = plt.subplots(figsize=(6, 5), dpi=160)
    distance_grid = transition["distance_grid"]
    if np.isfinite(distance_grid).any():
        vmax = float(np.nanquantile(distance_grid, 0.95))
        vmax = max(vmax, 1.0)
        image = ax.imshow(distance_grid, cmap="viridis", interpolation="nearest", vmin=0.0, vmax=vmax)
        colorbar_label = "Spot distance"
    else:
        background = fields["C_smooth_grid"] - fields["B_smooth_grid"]
        image = ax.imshow(background, cmap="coolwarm", interpolation="nearest")
        colorbar_label = "R = C - S"
        ax.text(
            0.5,
            0.5,
            "No interface detected\nTransition zone not defined",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=11,
            bbox={"facecolor": "white", "edgecolor": "0.4", "alpha": 0.85},
        )
    ty, tx = np.where(transition["transition_zone_grid"])
    if len(tx) > 0:
        ax.scatter(tx, ty, s=10, c="gold", alpha=0.65, marker="s", label="Transition zone")
    iy, ix = np.where(interface["interface_grid"])
    if len(ix) > 0:
        ax.scatter(ix, iy, s=10, c="black", marker="s", label="Interface")
    ax.set_title("SpatialTX transition-zone map")
    ax.set_axis_off()
    fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04, label=colorbar_label)
    if len(tx) > 0 or len(ix) > 0:
        ax.legend(loc="lower right", fontsize=7)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
