from __future__ import annotations

from dataclasses import asdict, is_dataclass
from typing import Any

GRAPH_SCHEMA_VERSION = "0.4"

FDR_SCOPE = "within_sample_within_analysis_table"
PERMUTATION_LIMITATION = (
    "Permutation P-values are exploratory and rely on exchangeability assumptions. "
    "They do not fully preserve the original spatial autocorrelation structure."
)
SMOOTHING_LIMITATION = (
    "Graph-smoothed context fields are intended for visualization and exploratory sensitivity analysis. "
    "Association statistics computed on fields smoothed over the same graph may be inflated and should "
    "not be interpreted as independent confirmatory evidence."
)
SPATIAL_ASSOCIATION_LIMITATION = (
    "These analyses provide exploratory spatial association and organization summaries. "
    "They do not establish causal, physical, or biological cell-cell interactions."
)

GRAPH_KEYS = {
    "radius": {
        "connectivities": "spatialtx_connectivities_radius",
        "distances": "spatialtx_distances_radius",
    },
    "knn": {
        "connectivities": "spatialtx_connectivities_knn",
        "distances": "spatialtx_distances_knn",
    },
    "lattice": {
        "connectivities": "spatialtx_connectivities_lattice",
        "distances": "spatialtx_distances_lattice",
    },
}


def json_safe(value: Any) -> Any:
    """Return a JSON-friendly representation for metadata dictionaries."""
    if is_dataclass(value):
        return json_safe(asdict(value))
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    try:
        import numpy as np

        if isinstance(value, np.generic):
            return value.item()
        if isinstance(value, np.ndarray):
            return value.tolist()
    except Exception:
        pass
    return value
