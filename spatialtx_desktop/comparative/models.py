from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from threading import Event
from typing import Callable

import pandas as pd


COMPARISON_MODES = ("pairwise", "paired", "unpaired", "manifest_batch")
NO_REGISTRATION_NOTICE = "Sample-level comparative summary; no direct spatial registration performed."
EXPLORATORY_NOTICE = (
    "Exploratory descriptive result only; not intended for diagnosis, treatment selection, "
    "or clinical decision-making."
)


@dataclass(frozen=True)
class SampleRecord:
    sample_id: str
    file_path: Path
    group: str
    pair_id: str = ""
    condition: str = ""
    batch: str = ""
    notes: str = ""

    def to_dict(self) -> dict:
        payload = asdict(self)
        payload["file_path"] = str(self.file_path)
        return payload


@dataclass
class ComparativeConfig:
    mode: str
    reference: str
    target: str
    c_genes: list[str]
    s_genes: list[str]
    c_q: float = 0.80
    s_q: float = 0.80
    g_q: float = 0.60
    scoring_options: dict = field(default_factory=dict)
    graph_settings: dict = field(default_factory=lambda: {"method": "knn", "k": 6, "weighting": "binary"})
    enable_h_expr: bool = True
    enable_v_expr: bool = True
    h_genes: list[str] | None = None
    v_genes: list[str] | None = None
    context_score_method: str = "z_score_mean"
    context_smoothing: str = "none"
    context_min_coverage: float = 0.25
    statistical_test: str = "auto"
    fdr_method: str = "benjamini-hochberg"
    seed: int = 42
    bootstrap_iterations: int = 1000
    use_cache: bool = True
    large_batch_warning_count: int = 50
    low_regime_confidence_threshold: float = 0.10

    def validate(self) -> None:
        if self.mode not in COMPARISON_MODES:
            raise ValueError(f"Unsupported comparison mode: {self.mode}")
        if not self.reference.strip() or not self.target.strip():
            raise ValueError("Reference and target definitions are required.")
        if self.reference.strip() == self.target.strip():
            raise ValueError("Reference and target must be different.")
        if not self.c_genes or not self.s_genes:
            raise ValueError("C-side and S-side gene programs are required.")
        for name, value in (("C", self.c_q), ("S", self.s_q), ("G", self.g_q)):
            if not 0 < float(value) < 1:
                raise ValueError(f"{name} quantile must be between 0 and 1.")
        if int(self.seed) < 0:
            raise ValueError("Random seed must be non-negative.")
        if int(self.bootstrap_iterations) < 100:
            raise ValueError("Bootstrap iterations must be at least 100.")
        if self.fdr_method != "benjamini-hochberg":
            raise ValueError("v0.6-beta supports Benjamini-Hochberg FDR only.")

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class ComparativeRunResult:
    run_dir: Path
    sample_metrics: pd.DataFrame
    delta_metrics: pd.DataFrame
    group_statistics: pd.DataFrame
    regime_transitions: pd.DataFrame
    warnings: pd.DataFrame
    run_manifest: pd.DataFrame
    figures: list[Path]
    summary_text: str
    effective_mode: str
    sample_scale: pd.DataFrame = field(default_factory=pd.DataFrame)
    metric_change_table: pd.DataFrame = field(default_factory=pd.DataFrame)
    normalized_metrics: pd.DataFrame = field(default_factory=pd.DataFrame)
    relative_changes: pd.DataFrame = field(default_factory=pd.DataFrame)
    scale_warnings: pd.DataFrame = field(default_factory=pd.DataFrame)
    hv_summary: pd.DataFrame = field(default_factory=pd.DataFrame)


Progress = Callable[[str], None]
CancelEvent = Event | None


class ComparativeCancelled(RuntimeError):
    """Raised after partial state is written when a comparative run is cancelled."""
