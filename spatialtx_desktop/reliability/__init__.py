"""Additive v0.65 Multi-axis Reliability Layer.

The reliability package preserves existing v0.6 signed C/S/R for legacy
Balance and consumes a separate pre-centering nonnegative program abundance
for Activity/Direction/Co-activation.  It never mutates FRAME2.6 C/S/R fields.
"""

from .core import compute_axis_reliability, compute_reliability_axes
from .models import (
    RELIABILITY_SCHEMA_VERSION,
    AxisReliabilityResult,
    ReliabilityConfig,
    ReliabilityResult,
)

__all__ = [
    "RELIABILITY_SCHEMA_VERSION",
    "AxisReliabilityResult",
    "ReliabilityConfig",
    "ReliabilityResult",
    "compute_axis_reliability",
    "compute_reliability_axes",
]
