"""SpatialTX v0.4 spatial graph and neighborhood analysis components.

Submodules are intentionally imported lazily. This keeps the existing desktop
application and Advanced Analysis UI importable even before optional sparse
scientific dependencies are checked, while the graph workflow itself still
requires the dependencies listed in ``requirements-desktop.txt``.
"""

__all__ = [
    "builder",
    "context",
    "continuous",
    "enrichment",
    "metadata",
    "plotting",
    "qc",
    "runner",
    "weights",
]
