from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Mapping

import numpy as np
import pandas as pd


RELIABILITY_SCHEMA_VERSION = "v0.65-reliability-v3-metric-qc"
RELIABILITY_STATES = (
    "low_activity",
    "c_dominant_active",
    "s_dominant_active",
    "active_coactivation_candidate",
)


@dataclass(frozen=True)
class ReliabilityConfig:
    """Configuration for the additive reliability sidecar.

    ``score_domain='nonnegative'`` is an interpretation requirement, not a
    transformation. Negative existing scores are retained in provenance and
    marked invalid for Activity/Direction/Co-activation calculations.
    """

    enabled: bool = False
    score_domain: str = "nonnegative"
    epsilon: float = 1.0e-9
    classification_enabled: bool = False
    activity_threshold: float | None = None
    direction_threshold: float | None = None
    strict_cross_exclusivity: bool = True
    dependence_qc: bool = True
    bootstrap_iterations: int = 1000
    permutation_iterations: int = 1000
    fdr_method: str = "benjamini-hochberg"
    seed: int = 42
    minimum_valid_spots: int = 30
    minimum_valid_fraction: float = 0.80
    warning_valid_fraction: float = 0.50
    dependence_abs_correlation_warning: float = 0.80
    gene_coverage_caution: float = 0.80
    gene_coverage_low: float = 0.50
    canonical_aliases: Mapping[str, str] = field(default_factory=dict)
    canonicalization_source: str = "user_supplied_or_none"
    canonicalization_version: str = "unversioned"

    def validate(self) -> None:
        if self.score_domain != "nonnegative":
            raise ValueError(
                "v0.65 Reliability Activity/Direction/Co-activation currently requires "
                "score_domain='nonnegative'."
            )
        if not np.isfinite(float(self.epsilon)) or float(self.epsilon) <= 0:
            raise ValueError("reliability epsilon must be a finite value greater than zero.")
        if int(self.bootstrap_iterations) < 1:
            raise ValueError("reliability bootstrap_iterations must be at least 1.")
        if int(self.permutation_iterations) < 1:
            raise ValueError("reliability permutation_iterations must be at least 1.")
        if int(self.seed) < 0:
            raise ValueError("reliability seed must be non-negative.")
        if int(self.minimum_valid_spots) < 3:
            raise ValueError("reliability minimum_valid_spots must be at least 3.")
        if not 0 <= float(self.warning_valid_fraction) <= float(self.minimum_valid_fraction) <= 1:
            raise ValueError(
                "reliability fractions must satisfy 0 <= warning_valid_fraction "
                "<= minimum_valid_fraction <= 1."
            )
        if self.fdr_method != "benjamini-hochberg":
            raise ValueError("v0.65 Reliability Layer supports Benjamini-Hochberg FDR only.")
        if not 0 < float(self.dependence_abs_correlation_warning) <= 1:
            raise ValueError("dependence_abs_correlation_warning must be in (0, 1].")
        if not 0 <= float(self.gene_coverage_low) <= float(self.gene_coverage_caution) <= 1:
            raise ValueError("gene coverage thresholds must satisfy 0 <= low <= caution <= 1.")
        if self.classification_enabled:
            if self.activity_threshold is None or self.direction_threshold is None:
                raise ValueError(
                    "classified mode requires explicit activity_threshold and direction_threshold."
                )
            if not np.isfinite(float(self.activity_threshold)) or float(self.activity_threshold) < 0:
                raise ValueError("activity_threshold must be a finite non-negative value.")
            if not np.isfinite(float(self.direction_threshold)) or not 0 <= float(self.direction_threshold) <= 1:
                raise ValueError("direction_threshold must be a finite value in [0, 1].")

    def to_dict(self) -> dict:
        payload = asdict(self)
        payload["canonical_aliases"] = dict(self.canonical_aliases)
        return payload

    @classmethod
    def from_value(cls, value: "ReliabilityConfig | Mapping | None") -> "ReliabilityConfig":
        if value is None:
            return cls()
        if isinstance(value, cls):
            return value
        if isinstance(value, Mapping):
            return cls(**dict(value))
        raise TypeError("reliability_config must be ReliabilityConfig, a mapping, or None.")


@dataclass(frozen=True)
class AxisReliabilityResult:
    axis: str
    # Existing v0.6 signed fields.  These remain the only source of legacy B.
    C: np.ndarray
    S: np.ndarray
    balance_B: np.ndarray
    balance_finite_input: np.ndarray
    # Separate pre-centering, nonnegative program abundance used by A/D/CA.
    activity_C: np.ndarray
    activity_S: np.ndarray
    activity_balance: np.ndarray
    activity_A: np.ndarray
    direction_D: np.ndarray
    ca_strength: np.ndarray
    ca_fraction: np.ndarray
    reliability_state: np.ndarray
    status: np.ndarray
    finite_input: np.ndarray
    nonnegative_input: np.ndarray
    valid_input: np.ndarray
    direction_defined: np.ndarray
    ca_defined: np.ndarray
    balance_score_source: str = "legacy_signed_cs"
    balance_score_domain: str = "signed"
    activity_score_source: str = "explicit_nonnegative_program_abundance"
    activity_score_domain: str = "nonnegative"
    activity_source_transformations: str = "explicit input; no clipping or shift"
    activity_source_version: str = "v0.65-nonnegative-program-mean-v1"


@dataclass
class ReliabilityResult:
    spot_results: pd.DataFrame = field(default_factory=pd.DataFrame)
    pair_summary: pd.DataFrame = field(default_factory=pd.DataFrame)
    metric_qc: pd.DataFrame = field(default_factory=pd.DataFrame)
    gene_coverage: pd.DataFrame = field(default_factory=pd.DataFrame)
    cross_exclusivity_audit: pd.DataFrame = field(default_factory=pd.DataFrame)
    axis_dependence_long: pd.DataFrame = field(default_factory=pd.DataFrame)
    axis_dependence_matrix: pd.DataFrame = field(default_factory=pd.DataFrame)
    score_domain_diagnostic: pd.DataFrame = field(default_factory=pd.DataFrame)
    score_domain_diagnostic_metadata: dict = field(default_factory=dict)
    qc: dict = field(default_factory=dict)
    files: list[Path] = field(default_factory=list)
