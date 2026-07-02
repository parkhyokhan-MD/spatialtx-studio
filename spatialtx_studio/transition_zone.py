from __future__ import annotations

import numpy as np
from scipy.ndimage import distance_transform_edt

from .spatial_fields import unlattice


def compute_transition_zone(fields: dict, interface: dict, config: dict) -> dict:
    interface_grid = interface["interface_grid"]
    tissue_mask = fields["tissue_mask_grid"]
    radius = float(config["transition_zone"]["radius_spots"])

    if interface_grid.sum() == 0:
        distance_grid = np.full(interface_grid.shape, np.nan, dtype=float)
        distance_spot = np.full(len(fields["rows"]), np.nan, dtype=float)
        transition_grid = np.zeros(interface_grid.shape, dtype=bool)
        transition_spot = np.zeros(len(fields["rows"]), dtype=bool)
    else:
        distance_grid = distance_transform_edt(~interface_grid)
        distance_grid = np.where(tissue_mask, distance_grid, np.nan)
        distance_spot = unlattice(np.nan_to_num(distance_grid, nan=-1.0), fields["rows"], fields["cols"])
        distance_spot[distance_spot < 0] = np.nan
        transition_grid = (distance_grid <= radius) & tissue_mask
        transition_spot = unlattice(transition_grid.astype(float), fields["rows"], fields["cols"]).astype(bool)

    return {
        "distance_grid": distance_grid,
        "distance_spot": distance_spot,
        "transition_zone_grid": transition_grid,
        "transition_zone_spot": transition_spot,
    }
