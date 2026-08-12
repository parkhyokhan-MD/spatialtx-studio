"""Comparative Spatial Transition Analysis for SpatialTX Studio v0.6-beta."""

from .models import ComparativeConfig, ComparativeRunResult, SampleRecord
from .multi_pair import (
    ComparabilityConfig,
    MultiPairRunResult,
    PairInterpretationConfig,
    PairSpec,
    run_multi_pair_analysis,
)
from .runner import run_comparative_analysis

__all__ = [
    "ComparativeConfig",
    "ComparativeRunResult",
    "ComparabilityConfig",
    "MultiPairRunResult",
    "PairInterpretationConfig",
    "PairSpec",
    "SampleRecord",
    "run_comparative_analysis",
    "run_multi_pair_analysis",
]
