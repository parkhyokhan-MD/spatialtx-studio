from __future__ import annotations

import os
import queue
import subprocess
import threading
import tkinter as tk
from dataclasses import asdict
from pathlib import Path
from threading import Event
from tkinter import filedialog, messagebox, ttk
from typing import Callable

import pandas as pd
from PIL import Image, ImageTk

from .comparative.models import ComparativeConfig
from .comparative.multi_pair import (
    MAX_MULTI_PAIRS,
    MULTI_PAIR_METRICS,
    ComparabilityConfig,
    MultiPairRunResult,
    PairInterpretationConfig,
    PairSpec,
    run_multi_pair_analysis,
    validate_pair_identity,
    validate_pair_specs,
)
from .reliability.models import ReliabilityConfig


def _reliability_summary_for_display(table: pd.DataFrame) -> pd.DataFrame:
    """Add human-readable metric support and N/A inference labels for the GUI."""

    display = table.copy()
    if display.empty:
        return display

    def is_true(value) -> bool:
        return value is True or str(value).strip().casefold() == "true"

    def support(row: pd.Series, metric: str) -> str:
        parts: list[str] = []
        for role in ("pre", "post"):
            defined = pd.to_numeric(
                pd.Series([row.get(f"{role}_{metric}_defined_n")]), errors="coerce"
            ).iloc[0]
            valid = pd.to_numeric(
                pd.Series([row.get(f"{role}_{metric}_valid_input_n")]), errors="coerce"
            ).iloc[0]
            fraction = pd.to_numeric(
                pd.Series([row.get(f"{role}_{metric}_defined_fraction")]), errors="coerce"
            ).iloc[0]
            if pd.notna(defined) and pd.notna(valid) and pd.notna(fraction):
                parts.append(
                    f"{role.title()} {int(defined)}/{int(valid)} ({float(fraction):.1%})"
                )
        return "; ".join(parts)

    def inference(row: pd.Series, metric: str, export_name: str) -> str:
        if not is_true(row.get(f"{metric}_inference_eligible", False)):
            reason = str(row.get(f"{metric}_qc_reason", "insufficient_metric_support"))
            return f"p: N/A — {reason}"
        p_value = pd.to_numeric(
            pd.Series([row.get(f"delta_{export_name}_permutation_p_value")]),
            errors="coerce",
        ).iloc[0]
        fdr = pd.to_numeric(
            pd.Series([row.get(f"delta_{export_name}_bh_fdr")]),
            errors="coerce",
        ).iloc[0]
        if pd.isna(p_value):
            return "p: N/A — eligible support but statistic undefined"
        fdr_text = "N/A" if pd.isna(fdr) else f"{float(fdr):.5g}"
        return f"p={float(p_value):.5g}; FDR={fdr_text}"

    display["direction_defined_spots"] = display.apply(
        lambda row: support(row, "direction"), axis=1
    )
    display["direction_inferential_test"] = display.apply(
        lambda row: inference(row, "direction", "D"), axis=1
    )
    display["ca_defined_spots"] = display.apply(lambda row: support(row, "ca"), axis=1)
    display["ca_inferential_test"] = display.apply(
        lambda row: inference(row, "ca", "CA_fraction"), axis=1
    )
    return display


class MultiPairAnalysisPanel(ttk.Frame):
    """One-to-six independent Pre/Post comparisons with a specimen comparability gate."""

    def __init__(
        self,
        parent,
        *,
        get_genes: Callable[[], tuple[list[str], list[str]]],
        get_quantiles: Callable[[], tuple[float, float, float]],
        get_scoring_options: Callable[[], object],
        get_output: Callable[[], str],
    ) -> None:
        super().__init__(parent, padding=8)
        self.get_genes = get_genes
        self.get_quantiles = get_quantiles
        self.get_scoring_options = get_scoring_options
        self.events: queue.Queue[tuple] = queue.Queue()
        self.cancel_event: Event | None = None
        self.busy = False
        self.last_result: MultiPairRunResult | None = None
        self.current_photo = None
        self._resize_job: str | None = None
        self.pair_vars: list[dict[str, tk.StringVar]] = []
        self.output_var = tk.StringVar(value=get_output())
        self.run_tag_var = tk.StringVar()
        self.graph_method_var = tk.StringVar(value="knn")
        self.graph_k_var = tk.IntVar(value=6)
        self.graph_radius_var = tk.StringVar()
        self.h_var = tk.BooleanVar(value=True)
        self.v_var = tk.BooleanVar(value=True)
        self.context_smoothing_var = tk.StringVar(value="none")
        self.reliability_layer_var = tk.BooleanVar(value=False)
        self.reliability_classification_var = tk.BooleanVar(value=False)
        self.reliability_epsilon_var = tk.StringVar(value="1e-9")
        self.reliability_activity_threshold_var = tk.StringVar()
        self.reliability_direction_threshold_var = tk.StringVar()
        self.status_var = tk.StringVar(
            value=f"Ready. Configure 1–{MAX_MULTI_PAIRS} complete Pre/Post pairs."
        )
        self.figure_var = tk.StringVar()
        self.figure_label_to_path: dict[str, Path] = {}
        self.comparability_config = ComparabilityConfig()
        self.interpretation_config = PairInterpretationConfig()
        self._build()
        self.after(100, self._poll)

    def _build(self) -> None:
        intro = ttk.LabelFrame(self, text="Paired Spatial Comparison — v0.65", padding=9)
        intro.pack(fill="x")
        ttk.Label(
            intro,
            text=(
                f"Runs up to {MAX_MULTI_PAIRS} independent Pre/Post pairs with shared C/S settings and optional H/V context axes. "
                "A transparent specimen "
                "comparability gate is reported separately from the observed spatial change. No registration or "
                "spot-wise subtraction is performed. Research use only; not for clinical decisions."
            ),
            foreground="#7c2d12",
            wraplength=1050,
            justify="left",
        ).pack(anchor="w", fill="x")

        input_frame = ttk.LabelFrame(self, text="1  Pair selection", padding=8)
        input_frame.pack(fill="x", pady=(7, 0))
        for column, text in enumerate(("Pair", "Label", "Site", "Pre H5AD", "", "Post H5AD", "", "")):
            ttk.Label(input_frame, text=text).grid(row=0, column=column, sticky="w", padx=(0, 5))
        for index in range(MAX_MULTI_PAIRS):
            variables = {
                "label": tk.StringVar(value=f"Pair_{index + 1}" if index == 0 else ""),
                "pre": tk.StringVar(),
                "post": tk.StringVar(),
                "site": tk.StringVar(value="unknown_site"),
            }
            self.pair_vars.append(variables)
            ttk.Label(input_frame, text=str(index + 1), width=4).grid(row=index + 1, column=0, sticky="w", pady=2)
            ttk.Entry(input_frame, textvariable=variables["label"], width=16).grid(row=index + 1, column=1, sticky="ew", padx=(0, 5), pady=2)
            ttk.Combobox(
                input_frame,
                textvariable=variables["site"],
                values=("same_site", "different_site", "unknown_site"),
                state="readonly",
                width=14,
            ).grid(row=index + 1, column=2, sticky="ew", padx=(0, 5), pady=2)
            pre_entry = ttk.Entry(input_frame, textvariable=variables["pre"], state="readonly")
            pre_entry.grid(row=index + 1, column=3, sticky="ew", padx=(0, 5), pady=2)
            ttk.Button(
                input_frame,
                text="Browse…",
                command=lambda i=index: self._browse_pair_file(i, "pre"),
            ).grid(row=index + 1, column=4, sticky="ew", padx=(0, 8), pady=2)
            post_entry = ttk.Entry(input_frame, textvariable=variables["post"], state="readonly")
            post_entry.grid(row=index + 1, column=5, sticky="ew", padx=(0, 5), pady=2)
            ttk.Button(
                input_frame,
                text="Browse…",
                command=lambda i=index: self._browse_pair_file(i, "post"),
            ).grid(row=index + 1, column=6, sticky="ew", pady=2)
            ttk.Button(
                input_frame,
                text="Clear",
                command=lambda i=index: self._clear_pair_row(i),
            ).grid(row=index + 1, column=7, sticky="ew", padx=(5, 0), pady=2)
            variables["pre_entry"] = pre_entry  # type: ignore[assignment]
            variables["post_entry"] = post_entry  # type: ignore[assignment]
        input_frame.columnconfigure(1, weight=1)
        input_frame.columnconfigure(3, weight=4)
        input_frame.columnconfigure(5, weight=4)

        settings = ttk.LabelFrame(self, text="2  Shared settings and output", padding=8)
        settings.pack(fill="x", pady=(7, 0))
        ttk.Label(
            settings,
            text="C/S programs, quantile thresholds, and scoring options are reused unchanged from Main Mapper for every pair.",
            foreground="#4b5563",
        ).grid(row=0, column=0, columnspan=8, sticky="w")
        ttk.Label(settings, text="Context graph").grid(row=1, column=0, sticky="w", pady=(5, 0))
        ttk.Combobox(
            settings,
            textvariable=self.graph_method_var,
            values=("knn", "radius", "lattice"),
            state="readonly",
            width=10,
        ).grid(row=2, column=0, sticky="ew")
        ttk.Label(settings, text="K").grid(row=1, column=1, sticky="w", padx=(6, 0), pady=(5, 0))
        ttk.Spinbox(settings, from_=1, to=50, textvariable=self.graph_k_var, width=6).grid(row=2, column=1, sticky="ew", padx=(6, 0))
        ttk.Label(settings, text="Radius (optional)").grid(row=1, column=2, sticky="w", padx=(6, 0), pady=(5, 0))
        ttk.Entry(settings, textvariable=self.graph_radius_var, width=10).grid(row=2, column=2, sticky="ew", padx=(6, 8))
        ttk.Label(settings, text="Run tag (optional)").grid(row=1, column=3, sticky="w", pady=(5, 0))
        ttk.Entry(settings, textvariable=self.run_tag_var, width=16).grid(row=2, column=3, sticky="ew", padx=(0, 8))
        ttk.Label(settings, text="Output root").grid(row=1, column=4, sticky="w", pady=(5, 0))
        ttk.Entry(settings, textvariable=self.output_var).grid(row=2, column=4, sticky="ew")
        ttk.Button(settings, text="Browse…", command=self._browse_output).grid(row=2, column=5, sticky="ew", padx=(5, 0))
        settings.columnconfigure(4, weight=1)
        context = ttk.Frame(settings)
        context.grid(row=3, column=0, columnspan=8, sticky="ew", pady=(7, 0))
        ttk.Checkbutton(
            context,
            text="H hypoxia-associated expression context",
            variable=self.h_var,
        ).pack(side="left")
        ttk.Checkbutton(
            context,
            text="V endothelial/angiogenic expression proxy",
            variable=self.v_var,
        ).pack(side="left", padx=(10, 0))
        ttk.Label(context, text="Context smoothing").pack(side="left", padx=(14, 4))
        ttk.Combobox(
            context,
            textvariable=self.context_smoothing_var,
            values=("none", "graph_mean"),
            state="readonly",
            width=11,
        ).pack(side="left")
        ttk.Label(
            context,
            text="Existing SpatialTX H/V programs; observational only.",
            foreground="#4b5563",
        ).pack(side="left", padx=(10, 0))
        reliability_controls = ttk.Frame(settings)
        reliability_controls.grid(row=4, column=0, columnspan=8, sticky="ew", pady=(7, 0))
        ttk.Checkbutton(
            reliability_controls,
            text="Enable v0.65 Reliability Layer (additive sidecars)",
            variable=self.reliability_layer_var,
        ).pack(side="left")
        ttk.Label(reliability_controls, text="epsilon").pack(side="left", padx=(12, 4))
        ttk.Entry(
            reliability_controls,
            textvariable=self.reliability_epsilon_var,
            width=9,
        ).pack(side="left")
        ttk.Checkbutton(
            reliability_controls,
            text="Classified mode",
            variable=self.reliability_classification_var,
        ).pack(side="left", padx=(12, 0))
        ttk.Label(reliability_controls, text="Activity threshold").pack(side="left", padx=(8, 4))
        ttk.Entry(
            reliability_controls,
            textvariable=self.reliability_activity_threshold_var,
            width=8,
        ).pack(side="left")
        ttk.Label(reliability_controls, text="Direction threshold").pack(side="left", padx=(8, 4))
        ttk.Entry(
            reliability_controls,
            textvariable=self.reliability_direction_threshold_var,
            width=8,
        ).pack(side="left")

        actions = ttk.Frame(self)
        actions.pack(fill="x", pady=(7, 0))
        self.run_button = ttk.Button(
            actions,
            text="Run Multi-Pair Comparison",
            style="Primary.TButton",
            command=self._run_async,
        )
        self.run_button.pack(side="left")
        self.cancel_button = ttk.Button(actions, text="Cancel safely", command=self._cancel, state="disabled")
        self.cancel_button.pack(side="left", padx=(6, 0))
        self.open_button = ttk.Button(actions, text="Open results folder", command=self._open_results, state="disabled")
        self.open_button.pack(side="left", padx=(6, 0))
        ttk.Button(actions, text="How to read results", command=self._show_rules_tab).pack(side="left", padx=(6, 0))
        self.progress = ttk.Progressbar(actions, mode="indeterminate", length=150)
        self.progress.pack(side="left", padx=(12, 6))
        ttk.Label(actions, textvariable=self.status_var, wraplength=720).pack(side="left", fill="x", expand=True)

        self.result_tabs = ttk.Notebook(self)
        self.result_tabs.pack(fill="both", expand=True, pady=(7, 0))
        interpretation_tab = ttk.Frame(self.result_tabs, padding=5)
        balance_tab = ttk.Frame(self.result_tabs, padding=5)
        spatial_tab = ttk.Frame(self.result_tabs, padding=5)
        context_tab = ttk.Frame(self.result_tabs, padding=5)
        reliability_tab = ttk.Frame(self.result_tabs, padding=5)
        axis_reliability_tab = ttk.Frame(self.result_tabs, padding=5)
        overview_tab = ttk.Frame(self.result_tabs, padding=5)
        figure_tab = ttk.Frame(self.result_tabs, padding=5)
        log_tab = ttk.Frame(self.result_tabs, padding=5)
        self.rules_tab = ttk.Frame(self.result_tabs, padding=7)
        self.result_tabs.add(interpretation_tab, text="Pair interpretation")
        self.result_tabs.add(balance_tab, text="1  Balance change")
        self.result_tabs.add(spatial_tab, text="2  Spatial organization")
        self.result_tabs.add(context_tab, text="H/V Context")
        self.result_tabs.add(reliability_tab, text="3  Specimen reliability")
        self.result_tabs.add(axis_reliability_tab, text="v0.65 Reliability")
        self.result_tabs.add(overview_tab, text="Multiaxial Overview")
        self.result_tabs.add(figure_tab, text="Figure")
        self.result_tabs.add(log_tab, text="Status log")
        self.result_tabs.add(self.rules_tab, text="Rules & interpretation")
        ttk.Label(
            interpretation_tab,
            text=(
                "Rule-based pair summary. Change magnitude and specimen reliability remain separate; this panel does "
                "not predict treatment response or combine the three layers into a clinical score."
            ),
            foreground="#7c2d12",
            wraplength=1100,
            justify="left",
        ).pack(anchor="w", fill="x")
        self.interpretation_tree = self._table(interpretation_tab, height=5)
        self.interpretation_tree.bind("<<TreeviewSelect>>", self._show_selected_pair_interpretation)
        self.interpretation_text = tk.Text(
            interpretation_tab,
            wrap="word",
            height=9,
            state="disabled",
            background="#f8fafc",
            padx=9,
            pady=7,
        )
        self.interpretation_text.pack(fill="both", expand=True, pady=(5, 0))
        ttk.Label(
            balance_tab,
            text="Layer 1 — Expression-derived program balance. C, S, and R=C-S are reported separately; this is not a spatial-organization or reliability score.",
            foreground="#1d4ed8",
            wraplength=1100,
            justify="left",
        ).pack(anchor="w", fill="x")
        self.balance_tree = self._table(balance_tab)
        self.pair_tree = self.balance_tree  # Backward-compatible UI attribute.
        ttk.Label(
            spatial_tab,
            text="Layer 2 — Coordinate-dependent spatial organization. Interface, diffuse, adjacency, fragmentation, and topology changes remain separate from C/S/R balance.",
            foreground="#047857",
            wraplength=1100,
            justify="left",
        ).pack(anchor="w", fill="x")
        self.spatial_tree = self._table(spatial_tab)
        ttk.Label(
            context_tab,
            text=(
                "Parallel observational context axes. H is hypoxia-associated expression context; V is an "
                "endothelial/angiogenic expression proxy, not perfusion or measured vascularity. "
                "Legacy within-sample centered-context q80 warnings are separate from the pair-pooled raw-context q90 summaries."
            ),
            foreground="#6b21a8",
            wraplength=1100,
            justify="left",
        ).pack(anchor="w", fill="x")
        self.context_tree = self._table(context_tab, height=6)
        self.context_tree.bind("<<TreeviewSelect>>", self._show_selected_context_summary)
        self.context_summary_text = tk.Text(
            context_tab,
            wrap="word",
            height=9,
            state="disabled",
            background="#faf5ff",
            padx=9,
            pady=7,
        )
        self.context_summary_text.pack(fill="x", pady=(5, 0))
        ttk.Label(
            reliability_tab,
            text="Layer 3 — Specimen reliability qualifies interpretation; it is not a biological outcome and is never combined with Layers 1 or 2 into one score.",
            foreground="#b45309",
            wraplength=1100,
            justify="left",
        ).pack(anchor="w", fill="x")
        self.reliability_tree = self._table(reliability_tab, height=4)
        ttk.Label(reliability_tab, text="Detailed technical, sampling, coverage, occupancy, and geometry checks:").pack(anchor="w", pady=(5, 0))
        self.qc_tree = self._table(reliability_tab, height=6)
        ttk.Label(
            axis_reliability_tab,
            text=(
                "Legacy signed Balance: preserved v0.6 C/S and B=C-S.  |  "
                "Nonnegative Activity/Co-activation: separate pre-centering program abundance; no clipping or shift. "
                "Validity requires both spot count and valid fraction; QC fail/warning reasons remain visible. "
                "Direction and CA_fraction require at least 30 defined spots and 80% of valid Activity inputs "
                "for inference; descriptive values remain visible when tests are not performed. "
                "Descriptive spot-distribution comparison of unregistered slides. Not specimen-level inference "
                "and not evidence of treatment effect."
            ),
            foreground="#0f766e",
            wraplength=1100,
            justify="left",
        ).pack(anchor="w", fill="x")
        reliability_notebook = ttk.Notebook(axis_reliability_tab)
        reliability_notebook.pack(fill="both", expand=True, pady=(5, 0))
        reliability_summary_tab = ttk.Frame(reliability_notebook, padding=4)
        reliability_audit_tab = ttk.Frame(reliability_notebook, padding=4)
        reliability_dependence_tab = ttk.Frame(reliability_notebook, padding=4)
        reliability_domain_tab = ttk.Frame(reliability_notebook, padding=4)
        reliability_notebook.add(reliability_summary_tab, text="Pre/Post summary")
        reliability_notebook.add(reliability_audit_tab, text="Gene audit & coverage")
        reliability_notebook.add(reliability_dependence_tab, text="Axis dependence")
        reliability_notebook.add(reliability_domain_tab, text="Score domain")
        self.axis_reliability_tree = self._table(reliability_summary_tab, height=10)
        self.reliability_audit_tree = self._table(reliability_audit_tab, height=10)
        self.axis_dependence_tree = self._table(reliability_dependence_tab, height=10)
        self.score_domain_tree = self._table(reliability_domain_tab, height=10)
        ttk.Label(
            overview_tab,
            text="The overview aligns the three result layers but does not calculate an overall response, quality, or clinical score.",
            foreground="#334155",
            wraplength=1100,
            justify="left",
        ).pack(anchor="w", fill="x")
        self.overview_tree = self._table(overview_tab, height=7)
        self.cohort_tree = self._table(overview_tab, height=7)
        figure_controls = ttk.Frame(figure_tab)
        figure_controls.pack(fill="x")
        self.figure_combo = ttk.Combobox(figure_controls, textvariable=self.figure_var, state="disabled")
        self.figure_combo.pack(side="left", fill="x", expand=True)
        self.figure_combo.bind("<<ComboboxSelected>>", lambda _event: self._show_figure())
        ttk.Button(figure_controls, text="Open figure", command=self._open_figure).pack(side="left", padx=(6, 0))
        self.figure_canvas = tk.Canvas(figure_tab, background="#111827", highlightthickness=0, height=270)
        self.figure_canvas.pack(fill="both", expand=True, pady=(5, 0))
        self.figure_canvas.bind("<Configure>", self._schedule_figure_resize)
        self.log_text = tk.Text(log_tab, wrap="word", state="disabled", background="#f8fafc")
        self.log_text.pack(fill="both", expand=True)
        self.rules_text = tk.Text(
            self.rules_tab,
            wrap="word",
            state="normal",
            background="#f8fafc",
            padx=10,
            pady=8,
        )
        self.rules_text.pack(fill="both", expand=True)
        self.rules_text.insert("1.0", self._rules_explanation())
        self.rules_text.configure(state="disabled")

    def _rules_explanation(self) -> str:
        cfg = self.comparability_config
        icfg = self.interpretation_config
        tolerances = {metric.export_name: metric.direction_tolerance for metric in MULTI_PAIR_METRICS}
        return (
            "HOW TO READ MULTI-PAIR RESULTS\n\n"
            "THREE SEPARATE RESULT LAYERS\n"
            "1. Balance change: C, S, and R=C-S.\n"
            "2. Spatial organization change: interface, diffuse, adjacency, fragmentation, transition burden, "
            "and compatible topology metrics.\n"
            "3. Specimen reliability: comparability plus technical, sampling, coverage, occupancy, and geometry QC.\n"
            "These layers are not combined into a single response, quality, or clinical score. Specimen reliability "
            "qualifies interpretation but is not itself a biological outcome or the probability that an interpretation is true.\n\n"
            "LAYER 1 — BALANCE CHANGE\n"
            "• Delta = Post - Pre. C, S, and R are shown separately.\n"
            "• Primary C/S/R state rows use existing field medians. Centered field means remain separate "
            "compatibility columns because they may be close to zero by construction.\n"
            "• Percent change is NA when the Pre value is zero/near-zero or when percentage interpretation is "
            "not meaningful for a signed centered score.\n"
            "• Direction: ↑ increase, ↓ decrease, → within the metric-specific near-zero tolerance.\n"
            f"  C/S/R tolerance={tolerances['C']:.3g}; fraction tolerance={tolerances['interface_fraction']:.3g}; "
            f"transition-burden tolerance={tolerances['transition_burden']:.3g}. Exact values remain in CSV.\n\n"
            "LAYER 2 — SPATIAL ORGANIZATION CHANGE\n"
            "• Spatial metrics require valid coordinates and describe organization rather than gene-program amount.\n"
            "• Interface, diffuse, adjacency, fragmentation, and topology changes remain separate measurements.\n"
            "• Type A/B/C and pattern transitions are derived exploratory descriptors. They are not specimen-reliability, "
            "treatment-response, or efficacy classifications.\n\n"
            "PARALLEL H/V CONTEXT AXES\n"
            "- H uses the existing hypoxia-associated expression-context program. V uses the existing "
            "endothelial/angiogenic expression proxy; V is not perfusion, vessel density, or measured vascularity.\n"
            "- H/V are optional. Missing H/V values are displayed as not available and never fail C/S/FRAME2.6 analysis.\n"
            "- H/V do not modify C, S, R, transition masks, interface/diffuse metrics, or Type A/B/C. They are not "
            "combined with the three result layers into a response score.\n\n"
            "OPTIONAL v0.65 RELIABILITY LAYER\n"
            "• Disabled by default. When disabled, v0.6 outputs and run metadata are unchanged.\n"
            "• Legacy signed Balance uses the preserved v0.6 C/S arrays. Activity/Direction/co-activation use a separate "
            "pre-z-score nonnegative program abundance from the same genes; neither source replaces the other.\n"
            "• Validity requires both the configured spot count and valid fraction. Negative, NaN, Inf, and zero-activity "
            "direction remain explicitly invalid or undefined; values remain visible with QC status.\n"
            "• Bootstrap/permutation outputs compare unregistered spot distributions only. They are not specimen-level "
            "inference or evidence of treatment effect.\n"
            "• Direction and CA_fraction inference requires metric-level PASS in both Pre and Post: at least 30 defined "
            "spots and at least 80% of valid Activity inputs. A 50–79.9% fraction is CAUTION; below 50% or fewer than "
            "30 spots is FAIL. Descriptive values remain visible, while unsupported CI/p/FDR remain N/A.\n"
            "• Continuous mode is primary. Classified mode requires explicit Activity and Direction thresholds and uses "
            "only low_activity, c_dominant_active, s_dominant_active, and active_coactivation_candidate.\n"
            "• Strict canonical gene overlap blocks the reliability run after writing its audit. Axis dependence is QC only; "
            "axes are never orthogonalized, transformed, or removed.\n"
            "• H/V are included in cross-exclusivity and coverage audit but are not treated as paired-pole reliability axes.\n\n"
            "TRANSPARENT QUALITATIVE CHANGE LABELS\n"
            "• Balance class uses max(|Delta C|, |Delta S|, |Delta R|): "
            f"Minimal < {icfg.balance_moderate_abs_delta:.3g}; Moderate < {icfg.balance_large_abs_delta:.3g}; Large otherwise.\n"
            "• Spatial class is Large when any configured metric reaches its Large threshold, Moderate when none is "
            "Large and at least one reaches its Moderate threshold, and Minimal otherwise.\n"
            f"  Interface Moderate/Large: {icfg.interface_moderate_abs_delta:.3g}/{icfg.interface_large_abs_delta:.3g}; "
            f"diffuse: {icfg.diffuse_moderate_abs_delta:.3g}/{icfg.diffuse_large_abs_delta:.3g}; "
            f"transition burden: {icfg.transition_burden_moderate_abs_delta:.3g}/{icfg.transition_burden_large_abs_delta:.3g}; "
            f"adjacency fractions: {icfg.adjacency_moderate_abs_delta:.3g}/{icfg.adjacency_large_abs_delta:.3g}; "
            f"fragmentation: {icfg.fragmentation_moderate_abs_delta:.3g}/{icfg.fragmentation_large_abs_delta:.3g}.\n"
            "• regime_preserved=yes means the broad candidate label is unchanged; it does not mean no biological change.\n"
            "• structure_preserved uses the regime result plus the spatial class. Low comparability downgrades a minimal, "
            "same-regime result from yes to probably.\n\n"
            "LAYER 3 — SPECIMEN RELIABILITY\n"
            "• Good: no configured primary warning and no secondary composition-proxy warning.\n"
            "• Caution: one primary caution, or a secondary composition-proxy warning.\n"
            f"• Low: any primary Low reason, or {cfg.caution_count_for_low} or more primary Caution reasons.\n"
            "• Missing optional QC is marked not_available and does not fail the batch.\n\n"
            "DEFAULT CAUTION / LOW THRESHOLDS\n"
            f"• Spot-count fold ratio: {cfg.spot_count_fold_caution:.2f} / {cfg.spot_count_fold_low:.2f}\n"
            f"• Feature-count fold ratio: {cfg.feature_count_fold_caution:.2f} / {cfg.feature_count_fold_low:.2f}\n"
            f"• Detected-genes fold ratio: {cfg.detected_genes_fold_caution:.2f} / {cfg.detected_genes_fold_low:.2f}\n"
            f"• Observed library-size fold ratio: {cfg.library_size_fold_caution:.2f} / {cfg.library_size_fold_low:.2f}\n"
            f"• In-tissue fraction absolute difference: {cfg.occupancy_difference_caution:.2f} / {cfg.occupancy_difference_low:.2f}\n"
            f"• Existing low-quality fraction difference: {cfg.low_quality_difference_caution:.2f} / {cfg.low_quality_difference_low:.2f}\n"
            f"• Valid-analysis-spot fraction difference: {cfg.valid_spot_difference_caution:.2f} / {cfg.valid_spot_difference_low:.2f}\n"
            f"• Spatial extent-area fold ratio: {cfg.extent_area_fold_caution:.2f} / {cfg.extent_area_fold_low:.2f}\n"
            f"• Tissue-component fold ratio: {cfg.tissue_component_fold_caution:.2f} / {cfg.tissue_component_fold_low:.2f}\n"
            f"• Required C/S gene coverage: below {cfg.gene_coverage_caution:.0%} / below {cfg.gene_coverage_low:.0%}\n"
            f"• Secondary C/S composition-proxy difference: {cfg.composition_proxy_difference_caution:.2f} / "
            f"{cfg.composition_proxy_difference_low:.2f}. This proxy cannot by itself produce Low.\n\n"
            "PAIR-ID SAFETY CHECK\n"
            "• Explicit patient/sample-like IDs in filenames are compared conservatively. A possible mismatch produces "
            "a visible warning asking the user to confirm the intended comparison; it does not automatically block the run "
            "or change the comparability class. Accessions alone are not treated as patient IDs.\n\n"
            "ANATOMICAL-SITE METADATA\n"
            "- same_site, different_site, and unknown_site are user-supplied comparison metadata.\n"
            "- different_site displays SITE-SHIFT WARNING because site-specific tissue composition or microenvironment "
            "may contribute to observed differences. It does not exclude the pair or redefine its FRAME2.6 result.\n"
            "- unknown_site remains visible as unknown and is not silently assumed to be the same site.\n\n"
            "INTERPRETATION LIMIT\n"
            "A Low result does not suppress numerical changes. It warns that technical quality, sampled tissue region, "
            "geometry, occupancy, or composition differences may substantially influence them. Regime transitions and "
            "direction counts are descriptive. The three layers must be reviewed side by side and must not be collapsed "
            "into a single overall score; they do not establish treatment response, therapeutic efficacy, responder "
            "status, drug sensitivity, or clinical benefit."
        )

    def _show_rules_tab(self) -> None:
        self.result_tabs.select(self.rules_tab)

    @staticmethod
    def _table(parent, height: int = 10) -> ttk.Treeview:
        frame = ttk.Frame(parent)
        frame.pack(fill="both", expand=True, pady=(4, 0))
        tree = ttk.Treeview(frame, show="headings", height=height)
        xscroll = ttk.Scrollbar(frame, orient="horizontal", command=tree.xview)
        yscroll = ttk.Scrollbar(frame, orient="vertical", command=tree.yview)
        tree.configure(xscrollcommand=xscroll.set, yscrollcommand=yscroll.set)
        tree.grid(row=0, column=0, sticky="nsew")
        yscroll.grid(row=0, column=1, sticky="ns")
        xscroll.grid(row=1, column=0, sticky="ew")
        frame.columnconfigure(0, weight=1)
        frame.rowconfigure(0, weight=1)
        return tree

    @staticmethod
    def _fill_tree(tree: ttk.Treeview, table: pd.DataFrame, preferred: list[str]) -> None:
        tree.delete(*tree.get_children())
        columns = [column for column in preferred if column in table.columns]
        if not columns:
            columns = list(table.columns[:16])
        tree.configure(columns=columns)
        for column in columns:
            tree.heading(column, text=column)
            tree.column(column, width=min(260, max(85, len(column) * 8)), stretch=False)
        for index, row in table.iterrows():
            values = []
            for column in columns:
                value = row[column]
                if pd.isna(value):
                    value = ""
                elif isinstance(value, float):
                    value = f"{value:.5g}"
                values.append(value)
            tree.insert("", "end", iid=f"row_{index}", values=values)

    def _browse_pair_file(self, index: int, slot: str) -> None:
        value = filedialog.askopenfilename(
            title=f"Pair {index + 1} {slot.title()} H5AD",
            filetypes=[("AnnData H5AD", "*.h5ad")],
            parent=self,
        )
        if not value:
            return
        self.pair_vars[index][slot].set(value)
        entry = self.pair_vars[index].get(f"{slot}_entry")
        if hasattr(entry, "xview_moveto"):
            self.after_idle(lambda widget=entry: widget.xview_moveto(1.0))
        if not self.pair_vars[index]["label"].get().strip():
            self.pair_vars[index]["label"].set(f"Pair_{index + 1}")

    def _browse_output(self) -> None:
        value = filedialog.askdirectory(title="Multi-Pair output root", parent=self)
        if value:
            self.output_var.set(value)

    def _clear_pair_row(self, index: int) -> None:
        self.pair_vars[index]["label"].set("")
        self.pair_vars[index]["pre"].set("")
        self.pair_vars[index]["post"].set("")
        self.pair_vars[index]["site"].set("unknown_site")

    def _pairs(self) -> list[PairSpec]:
        pairs: list[PairSpec] = []
        partial: list[str] = []
        for index, variables in enumerate(self.pair_vars, start=1):
            label = variables["label"].get().strip()
            pre = variables["pre"].get().strip()
            post = variables["post"].get().strip()
            site = variables["site"].get().strip() or "unknown_site"
            if not pre and not post:
                if label and label != f"Pair_{index}":
                    partial.append(f"Pair {index} has a label but no Pre/Post files")
                continue
            if not pre or not post:
                partial.append(f"Pair {index} is partially filled; both Pre and Post are required")
                continue
            pairs.append(PairSpec(label or f"Pair_{index}", Path(pre), Path(post), site))
        if partial:
            raise ValueError("\n".join(partial))
        return validate_pair_specs(pairs)

    def _analysis_config(self) -> ComparativeConfig:
        c_genes, s_genes = self.get_genes()
        c_q, s_q, g_q = self.get_quantiles()
        scoring = self.get_scoring_options()
        scoring_options = asdict(scoring) if hasattr(scoring, "__dataclass_fields__") else dict(scoring)
        graph = {"method": self.graph_method_var.get(), "k": int(self.graph_k_var.get()), "weighting": "binary"}
        radius = self.graph_radius_var.get().strip()
        if radius:
            graph["radius"] = float(radius)
        return ComparativeConfig(
            mode="pairwise",
            reference="Pre",
            target="Post",
            c_genes=c_genes,
            s_genes=s_genes,
            c_q=c_q,
            s_q=s_q,
            g_q=g_q,
            scoring_options=scoring_options,
            graph_settings=graph,
            enable_h_expr=bool(self.h_var.get()),
            enable_v_expr=bool(self.v_var.get()),
            context_smoothing=self.context_smoothing_var.get(),
        )

    def _reliability_config(self) -> ReliabilityConfig:
        if not self.reliability_layer_var.get():
            return ReliabilityConfig(enabled=False)
        classified = bool(self.reliability_classification_var.get())
        activity_text = self.reliability_activity_threshold_var.get().strip()
        direction_text = self.reliability_direction_threshold_var.get().strip()
        config = ReliabilityConfig(
            enabled=True,
            score_domain="nonnegative",
            epsilon=float(self.reliability_epsilon_var.get().strip()),
            classification_enabled=classified,
            activity_threshold=float(activity_text) if activity_text else None,
            direction_threshold=float(direction_text) if direction_text else None,
            strict_cross_exclusivity=True,
            dependence_qc=True,
            bootstrap_iterations=1000,
            permutation_iterations=1000,
            fdr_method="benjamini-hochberg",
            seed=42,
        )
        config.validate()
        return config

    def _run_async(self) -> None:
        if self.busy:
            return
        try:
            pairs = self._pairs()
            config = self._analysis_config()
            reliability_config = self._reliability_config()
            output = self.output_var.get().strip()
            if not output:
                raise ValueError("Choose an output root.")
            config.validate()
        except Exception as exc:
            messagebox.showwarning("Multi-Pair Pre/Post", str(exc), parent=self)
            return
        pair_id_warnings = [
            validate_pair_identity(pair.pre_path, pair.post_path)["pair_id_warning"]
            for pair in pairs
        ]
        pair_id_warnings = [warning for warning in pair_id_warnings if warning]
        if pair_id_warnings and not messagebox.askokcancel(
            "Possible pair-ID mismatch",
            "\n\n".join(pair_id_warnings)
            + "\n\nPlease confirm these are intended comparisons. Continue with the run?",
            parent=self,
        ):
            self.status_var.set("Run cancelled so pair assignments can be reviewed.")
            return
        self.busy = True
        self.cancel_event = Event()
        self.progress.start(12)
        self.run_button.configure(state="disabled")
        self.cancel_button.configure(state="normal")
        self.status_var.set(f"Running {len(pairs)} independent Pre/Post pair(s)…")

        def progress(message: str) -> None:
            self.events.put(("progress", message))

        def work() -> None:
            try:
                result = run_multi_pair_analysis(
                    pairs,
                    config,
                    output,
                    self.comparability_config,
                    progress,
                    self.cancel_event,
                    self.run_tag_var.get().strip() or None,
                    interpretation_config=self.interpretation_config,
                    reliability_config=reliability_config,
                )
                self.events.put(("done", result))
            except Exception as exc:
                self.events.put(("error", exc))

        threading.Thread(target=work, daemon=True).start()

    def _cancel(self) -> None:
        if self.cancel_event is not None:
            self.cancel_event.set()
            self.status_var.set("Cancellation requested; the current specimen will finish before partial export.")

    def _finish_busy(self) -> None:
        self.busy = False
        self.cancel_event = None
        self.progress.stop()
        self.run_button.configure(state="normal")
        self.cancel_button.configure(state="disabled")

    def _poll(self) -> None:
        try:
            while True:
                event = self.events.get_nowait()
                if event[0] == "progress":
                    self.status_var.set(str(event[1]))
                elif event[0] == "done":
                    self._finish_busy()
                    self._load_result(event[1])
                elif event[0] == "error":
                    self._finish_busy()
                    self.status_var.set(str(event[1]))
                    messagebox.showerror("Multi-Pair Pre/Post failed", str(event[1]), parent=self)
        except queue.Empty:
            pass
        if self.winfo_exists():
            self.after(100, self._poll)

    def _load_result(self, result: MultiPairRunResult) -> None:
        self.last_result = result
        self._fill_tree(
            self.interpretation_tree,
            result.pair_interpretation_summary,
            [
                "pair_label", "balance_change_class", "spatial_change_class", "regime_transition",
                "regime_preserved", "structure_preserved", "comparability", "pair_status_flag",
                "interpretation_confidence", "site_comparability", "interpretive_flag",
                "pair_id_validation", "pair_id_warning", "status", "error",
            ],
        )
        self._fill_tree(
            self.balance_tree,
            result.balance_changes,
            [
                "pair_label", "status", "pre_C", "post_C", "delta_C", "direction_C",
                "pre_S", "post_S", "delta_S", "direction_S", "pre_R", "post_R", "delta_R", "direction_R",
                "percent_change_C", "percent_change_status_C", "percent_change_S", "percent_change_status_S",
                "percent_change_R", "percent_change_status_R", "error",
            ],
        )
        self._fill_tree(
            self.spatial_tree,
            result.spatial_organization_changes,
            [
                "pair_label", "status",
                "pre_interface_fraction", "post_interface_fraction", "delta_interface_fraction",
                "direction_interface_fraction",
                "pre_diffuse_fraction", "post_diffuse_fraction", "delta_diffuse_fraction",
                "direction_diffuse_fraction",
                "pre_transition_burden", "post_transition_burden", "delta_transition_burden",
                "direction_transition_burden",
                "delta_adj_same_fraction", "direction_adj_same_fraction",
                "delta_adj_zero_fraction", "direction_adj_zero_fraction",
                "delta_adj_opposite_fraction", "direction_adj_opposite_fraction",
                "delta_interface_fragmentation", "direction_interface_fragmentation",
                "regime_transition", "pattern_transition", "error",
            ],
        )
        self._fill_tree(
            self.context_tree,
            result.context_changes,
            [
                "pair_label", "status",
                "pre_H", "post_H", "delta_H", "pre_H_q90", "post_H_q90",
                "pre_H_high_fraction", "post_H_high_fraction",
                "pre_H_local_hotspot_fraction", "post_H_local_hotspot_fraction",
                "pre_H_context_status", "post_H_context_status",
                "pre_V", "post_V", "delta_V", "pre_V_q90", "post_V_q90",
                "pre_V_high_fraction", "post_V_high_fraction",
                "pre_V_local_hotspot_fraction", "post_V_local_hotspot_fraction",
                "pre_V_context_status", "post_V_context_status", "error",
            ],
        )
        self._fill_tree(
            self.reliability_tree,
            result.comparative_qc_summary,
            [
                "pair_id", "comparability_status", "interpretation_confidence", "site_comparability",
                "primary_mismatch_1", "primary_mismatch_2", "primary_mismatch_3",
                "comparability_reason", "site_shift_warning",
            ],
        )
        self._fill_tree(
            self.qc_tree,
            result.comparability_qc,
            [
                "pair_label", "comparability", "qc_metric", "category", "pre_value", "post_value",
                "comparison_value", "availability_status", "severity", "primary_for_classification", "reason",
            ],
        )
        reliability_summary_display = _reliability_summary_for_display(
            result.reliability_pair_summary
        )
        self._fill_tree(
            self.axis_reliability_tree,
            reliability_summary_display,
            [
                "pair_label", "axis", "balance_score_source", "balance_score_domain",
                "activity_score_source", "activity_score_domain",
                "pre_B", "post_B", "delta_B",
                "pre_A", "post_A", "delta_A", "pre_D", "post_D", "delta_D",
                "direction_defined_spots", "direction_qc_status",
                "direction_inference_eligible", "direction_inferential_test",
                "pre_CA_strength", "post_CA_strength", "delta_CA_strength",
                "pre_CA_fraction", "post_CA_fraction", "delta_CA_fraction",
                "ca_defined_spots", "ca_qc_status",
                "ca_inference_eligible", "ca_inferential_test",
                "pre_valid_input_fraction", "post_valid_input_fraction",
                "pre_undefined_fraction", "post_undefined_fraction",
                "pre_total_spot_count", "post_total_spot_count",
                "pre_valid_spot_count", "post_valid_spot_count", "minimum_valid_spots",
                "minimum_valid_fraction", "warning_valid_fraction",
                "pre_score_validity", "post_score_validity", "pair_score_validity",
                "pre_validity_reason", "post_validity_reason", "pair_validity_reason",
                "activity_summary_included_in_conclusion", "classification_enabled", "epsilon",
                "inference_level", "specimen_level_inference", "inference_warning",
            ],
        )
        exclusivity_display = result.cross_exclusivity_audit.copy()
        if not exclusivity_display.empty:
            exclusivity_display.insert(0, "audit_type", "cross_exclusivity")
        coverage_display = result.reliability_gene_coverage.copy()
        if not coverage_display.empty:
            coverage_display.insert(0, "audit_type", "gene_coverage")
        reliability_audit = pd.concat(
            [table for table in (exclusivity_display, coverage_display) if not table.empty],
            ignore_index=True,
            sort=False,
        ) if not exclusivity_display.empty or not coverage_display.empty else pd.DataFrame()
        self._fill_tree(
            self.reliability_audit_tree,
            reliability_audit,
            [
                "audit_type", "pair_label", "sample_role", "axis", "pole",
                "canonical_gene", "overlap_type", "severity", "action",
                "n_genes_requested", "n_genes_present", "gene_coverage_fraction",
                "missing_genes", "coverage_status", "score_validity", "normalization_rule",
            ],
        )
        self._fill_tree(
            self.axis_dependence_tree,
            result.axis_dependence_long,
            [
                "sample_id", "axis_i", "axis_j", "dependence_type", "metric_i", "metric_j",
                "pearson_correlation", "spearman_correlation", "valid_spot_count",
                "missing_undefined_fraction", "permutation_p_value", "bh_fdr", "qc_status",
                "metric_inference_eligible", "metric_inference_qc_reason",
                "permutation_scope",
            ],
        )
        self._fill_tree(
            self.score_domain_tree,
            result.reliability_score_domain_diagnostic,
            [
                "pair_label", "sample_role", "score_role", "score_source", "score_domain",
                "total_spots", "finite_spots", "C_min", "C_q01", "C_median", "C_q99", "C_max",
                "S_min", "S_q01", "S_median", "S_q99", "S_max",
                "either_negative_count", "either_negative_fraction",
                "both_nonnegative_count", "both_nonnegative_fraction",
                "nonfinite_count", "nonfinite_fraction", "first_negative_stage",
                "transformation_history",
            ],
        )
        self._fill_tree(
            self.overview_tree,
            result.overview_interpretation,
            [
                "pair_label", "comparability", "balance_change_class", "spatial_change_class",
                "regime_transition", "regime_preserved", "structure_preserved", "delta_C", "delta_S",
                "delta_R", "delta_interface", "delta_diffuse", "delta_burden", "delta_H", "delta_V",
                "interpretation_confidence", "site_comparability", "interpretive_flag",
                "pair_status_flag", "pair_id_warning", "site_shift_warning",
            ],
        )
        self._fill_tree(self.cohort_tree, result.cohort_summary, ["summary_metric", "count", "eligible_pairs"])
        self.log_text.configure(state="normal")
        self.log_text.delete("1.0", "end")
        self.log_text.insert("1.0", result.summary_text + "\n\n")
        for _, row in result.pair_interpretation_summary.iterrows():
            self.log_text.insert(
                "end",
                f"{row['pair_label']}: {row['status']} | Balance {row['balance_change_class']}"
                f" | Spatial {row['spatial_change_class']} | Reliability {row['comparability']}"
                f" | {row['interpretive_flag']}"
                f"{(' | ' + str(row['pair_id_warning'])) if str(row.get('pair_id_warning', '')).strip() else ''}"
                f"{(' | ERROR: ' + str(row['error'])) if str(row.get('error', '')).strip() else ''}\n",
            )
        self.log_text.configure(state="disabled")
        interpretation_rows = self.interpretation_tree.get_children()
        if interpretation_rows:
            self.interpretation_tree.selection_set(interpretation_rows[0])
            self.interpretation_tree.focus(interpretation_rows[0])
            self._show_selected_pair_interpretation()
        context_rows = self.context_tree.get_children()
        if context_rows:
            self.context_tree.selection_set(context_rows[0])
            self.context_tree.focus(context_rows[0])
            self._show_selected_context_summary()
        self.figure_label_to_path = {
            (
                "Three-layer overview"
                if path.name == "multi_pair_comparative_overview.png"
                else "Multiaxial change profile"
                if path.name == "multiaxial_pair_overview.png"
                else "v0.65 axis-dependence QC"
                if path.name == "axis_dependence_heatmap.png"
                else path.stem
            ): path
            for path in result.figures
        }
        labels = list(self.figure_label_to_path)
        self.figure_combo.configure(values=labels, state="readonly" if labels else "disabled")
        if labels:
            self.figure_var.set(labels[0])
            self._show_figure()
        self.open_button.configure(state="normal")
        passed = int(result.pair_results["status"].eq("PASS").sum())
        self.status_var.set(f"Multi-Pair complete: {passed}/{len(result.pair_results)} pairs passed. Results: {result.run_dir}")

    def _show_selected_pair_interpretation(self, _event=None) -> None:
        self.interpretation_text.configure(state="normal")
        self.interpretation_text.delete("1.0", "end")
        if self.last_result is None:
            self.interpretation_text.configure(state="disabled")
            return
        selection = self.interpretation_tree.selection()
        if not selection:
            self.interpretation_text.configure(state="disabled")
            return
        pair_label = self.interpretation_tree.set(selection[0], "pair_label")
        selected = self.last_result.pair_interpretation_summary.loc[
            self.last_result.pair_interpretation_summary["pair_label"].astype(str).eq(pair_label)
        ]
        if selected.empty:
            self.interpretation_text.configure(state="disabled")
            return
        row = selected.iloc[0]
        text = str(row.get("interpretation_text", ""))
        qc_context = [
            ("Technical mismatch", row.get("technical_mismatch_reason", "")),
            ("Sampling mismatch", row.get("sampling_mismatch_reason", "")),
            ("Composition proxy", row.get("composition_proxy_reason", "")),
            ("Pair-ID validation", row.get("pair_id_warning", "")),
            ("Site", row.get("site_shift_warning", "")),
        ]
        qc_summary = self.last_result.comparative_qc_summary
        if not qc_summary.empty:
            selected_qc = qc_summary.loc[qc_summary["pair_id"].astype(str).eq(pair_label)]
            if not selected_qc.empty:
                qc_context.append(("Primary mismatches", selected_qc.iloc[0].get("primary_mismatch_summary", "")))
        visible_context = [f"{label}: {value}" for label, value in qc_context if str(value).strip()]
        if visible_context:
            text += "\n\nQC context:\n" + "\n".join(visible_context)
        self.interpretation_text.insert("1.0", text)
        self.interpretation_text.configure(state="disabled")

    def _show_selected_context_summary(self, _event=None) -> None:
        self.context_summary_text.configure(state="normal")
        self.context_summary_text.delete("1.0", "end")
        if self.last_result is None:
            self.context_summary_text.configure(state="disabled")
            return
        selection = self.context_tree.selection()
        if not selection:
            self.context_summary_text.configure(state="disabled")
            return
        pair_label = self.context_tree.set(selection[0], "pair_label")
        selected = self.last_result.multiaxial_pair_summary.loc[
            self.last_result.multiaxial_pair_summary["pair_id"].astype(str).eq(pair_label)
        ]
        if selected.empty:
            self.context_summary_text.configure(state="disabled")
            return
        row = selected.iloc[0]

        def number(value) -> str:
            numeric = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
            return f"{float(numeric):.4g}" if pd.notna(numeric) else "NA"

        def fraction(value) -> str:
            numeric = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
            return f"{float(numeric):.1%}" if pd.notna(numeric) else "NA"

        lines = [f"Pair: {pair_label}"]
        for axis in ("H", "V"):
            lines.extend([
                "",
                f"{axis} observational context",
                f"  status: {row.get(f'{axis}_context_status_pre', 'NA')} -> {row.get(f'{axis}_context_status_post', 'NA')}",
                f"  gene coverage: {fraction(row.get(f'{axis}_gene_coverage_pre'))} -> {fraction(row.get(f'{axis}_gene_coverage_post'))}",
                f"  median: {number(row.get(f'{axis}_pre'))} -> {number(row.get(f'{axis}_post'))}",
                f"  q90: {number(row.get(f'{axis}_q90_pre'))} -> {number(row.get(f'{axis}_q90_post'))}",
                f"  pair-pooled high fraction: {fraction(row.get(f'{axis}_high_fraction_pre'))} -> {fraction(row.get(f'{axis}_high_fraction_post'))}",
                f"  high-context local fraction: {fraction(row.get(f'{axis}_local_fraction_pre'))} -> {fraction(row.get(f'{axis}_local_fraction_post'))}",
                f"  shared within-pair q90 threshold: {number(row.get(f'{axis}_pair_pooled_q90'))}",
            ])
        lines.extend([
            "",
            "Median is one summary only. A stable zero median does not erase a change in upper-tail or local high-context fractions.",
            "Audit note: a legacy single-sample centered-context q80 warning is not a warning about the pair-pooled raw-context q90 high fraction shown above.",
            "H/V remain observational and do not modify C/S/R, transition masks, comparability QC, or Type A/B/C.",
        ])
        self.context_summary_text.insert("1.0", "\n".join(lines))
        self.context_summary_text.configure(state="disabled")

    def _current_figure(self) -> Path | None:
        return self.figure_label_to_path.get(self.figure_var.get())

    def _schedule_figure_resize(self, _event=None) -> None:
        if self._resize_job is not None:
            try:
                self.after_cancel(self._resize_job)
            except tk.TclError:
                pass
        self._resize_job = self.after(140, self._show_figure)

    def _show_figure(self) -> None:
        self._resize_job = None
        self.figure_canvas.delete("all")
        path = self._current_figure()
        if path is None or not path.is_file():
            self.figure_canvas.create_text(20, 20, anchor="nw", fill="#cbd5e1", text="Multi-Pair figure unavailable.")
            self.current_photo = None
            return
        try:
            with Image.open(path) as source:
                display = source.convert("RGBA")
            display.thumbnail(
                (max(100, self.figure_canvas.winfo_width() - 20), max(100, self.figure_canvas.winfo_height() - 20)),
                Image.Resampling.LANCZOS,
            )
            self.current_photo = ImageTk.PhotoImage(display, master=self.figure_canvas)
            self.figure_canvas.create_image(
                max(1, self.figure_canvas.winfo_width()) // 2,
                max(1, self.figure_canvas.winfo_height()) // 2,
                image=self.current_photo,
            )
        except (OSError, ValueError) as exc:
            self.current_photo = None
            self.figure_canvas.create_text(20, 20, anchor="nw", fill="#cbd5e1", text=f"Could not display figure: {exc}")

    @staticmethod
    def _open_path(path: Path) -> None:
        if os.name == "nt":
            os.startfile(path)  # type: ignore[attr-defined]
        else:
            subprocess.Popen(["open" if os.uname().sysname == "Darwin" else "xdg-open", str(path)])

    def _open_results(self) -> None:
        if self.last_result and self.last_result.run_dir.is_dir():
            try:
                self._open_path(self.last_result.run_dir)
            except OSError as exc:
                messagebox.showerror("Open Multi-Pair results", str(exc), parent=self)

    def _open_figure(self) -> None:
        path = self._current_figure()
        if path and path.is_file():
            try:
                self._open_path(path)
            except OSError as exc:
                messagebox.showerror("Open Multi-Pair figure", str(exc), parent=self)
