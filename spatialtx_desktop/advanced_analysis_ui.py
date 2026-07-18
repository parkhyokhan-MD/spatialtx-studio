from __future__ import annotations

import os
import queue
import subprocess
import threading
import tkinter as tk
from pathlib import Path
from threading import Event
from tkinter import filedialog, messagebox, ttk
from typing import Callable

import numpy as np
import pandas as pd

from .advanced_analysis import MODULE_LABELS, run_advanced_batch
from .graph.context import DEFAULT_HYPOXIA_GENES, DEFAULT_VASCULAR_PROXY_GENES
from .workflow import parse_gene_text


class AdvancedAnalysisPanel(ttk.Frame):
    def __init__(
        self,
        parent,
        get_samples: Callable[[], list[Path]],
        get_genes: Callable[[], tuple[list[str], list[str]]],
        get_quantiles: Callable[[], tuple[float, float, float]],
        get_output: Callable[[], str],
    ) -> None:
        super().__init__(parent)
        self.get_samples = get_samples
        self.get_genes = get_genes
        self.get_quantiles = get_quantiles
        self.output_var = tk.StringVar(value=get_output())
        self.permutations_var = tk.IntVar(value=499)
        self.seed_var = tk.IntVar(value=20260705)
        self.graph_method_var = tk.StringVar(value="radius")
        self.graph_radius_var = tk.StringVar(value="")
        self.graph_coordinate_unit_var = tk.StringVar(value="native")
        self.graph_coordinate_scale_var = tk.StringVar(value="")
        self.graph_scale_source_var = tk.StringVar(value="")
        self.graph_k_var = tk.IntVar(value=6)
        self.graph_weighting_var = tk.StringVar(value="binary")
        self.graph_sym_var = tk.StringVar(value="union")
        self.graph_enable_h_var = tk.BooleanVar(value=True)
        self.graph_enable_v_var = tk.BooleanVar(value=True)
        self.graph_h_method_var = tk.StringVar(value="z_score_mean")
        self.graph_v_method_var = tk.StringVar(value="z_score_mean")
        self.graph_smoothing_var = tk.StringVar(value="none")
        self.graph_h_quantile_var = tk.StringVar(value="0.80")
        self.graph_v_quantile_var = tk.StringVar(value="0.80")
        self.graph_permutations_var = tk.IntVar(value=999)
        self.graph_seed_var = tk.IntVar(value=20260713)
        self.graph_min_coverage_var = tk.StringVar(value="0.25")
        self.graph_min_spot_fraction_var = tk.StringVar(value="0.01")
        self.graph_allow_low_coverage_var = tk.BooleanVar(value=False)
        self.graph_h_genes_var = tk.StringVar(value=", ".join(DEFAULT_HYPOXIA_GENES))
        self.graph_v_genes_var = tk.StringVar(value=", ".join(DEFAULT_VASCULAR_PROXY_GENES))
        self.graph_label_source_var = tk.StringVar(value="")
        self.graph_label_mode_var = tk.StringVar(value="categorical_state")
        self.graph_mask_a_var = tk.StringVar(value="")
        self.graph_mask_b_var = tk.StringVar(value="")
        self.graph_x_column_var = tk.StringVar(value="")
        self.graph_y_column_var = tk.StringVar(value="")
        self.graph_x_mode_var = tk.StringVar(value="continuous_score")
        self.graph_y_mode_var = tk.StringVar(value="continuous_score")
        self.graph_permutation_scope_var = tk.StringVar(value="whole_slide")
        self.graph_strata_column_var = tk.StringVar(value="")
        self.graph_tissue_only_var = tk.BooleanVar(value=False)
        self.graph_robustness_var = tk.BooleanVar(value=False)
        self.graph_cancel_event: Event | None = None
        self.events: queue.Queue[tuple] = queue.Queue()
        self.busy = False
        self.last_run: Path | None = None
        self.dashboard_records: dict[str, dict] = {}
        self.selected_record: dict | None = None
        self.dashboard_counter = 0
        self._build()
        self.graph_coordinate_unit_var.trace_add("write", self._update_coordinate_status)
        self.graph_coordinate_scale_var.trace_add("write", self._update_coordinate_status)
        self.graph_scale_source_var.trace_add("write", self._update_coordinate_status)
        self.after(100, self._poll)

    def _build(self) -> None:
        intro = ttk.LabelFrame(self, text="Advanced Analysis", padding=12)
        intro.pack(fill="x")
        ttk.Label(
            intro,
            text=(
                "Quantitative extensions of the unchanged Cx/Sx framework. Each module uses the samples, gene programs, "
                "and interface thresholds currently shown in the main workspace. Existing Transition Mapper outputs are not changed."
            ),
            wraplength=900,
            justify="left",
        ).pack(anchor="w", fill="x")
        output = ttk.Frame(intro)
        output.pack(fill="x", pady=(9, 0))
        ttk.Label(output, text="Output root").pack(side="left")
        ttk.Entry(output, textvariable=self.output_var).pack(side="left", fill="x", expand=True, padx=6)
        ttk.Button(output, text="Browse...", command=self._browse_output).pack(side="left")
        ttk.Button(output, text="Use main output", command=lambda: self.output_var.set(self.get_output())).pack(side="left", padx=(5, 0))
        run_all = ttk.Frame(intro)
        run_all.pack(fill="x", pady=(9, 0))
        ttk.Button(
            run_all,
            text="Run All 3 Analyses + Show Dashboard",
            style="Primary.TButton",
            command=lambda: self._run("all"),
        ).pack(side="left", fill="x", expand=True)
        ttk.Label(
            run_all,
            text="Runs Composition → Enrichment → Interaction for every selected sample.",
            foreground="#4b5563",
        ).pack(side="left", padx=(10, 0))

        self.tabs = ttk.Notebook(self)
        self.tabs.pack(fill="both", expand=True, pady=(10, 0))
        composition = ttk.Frame(self.tabs, padding=12)
        enrichment = ttk.Frame(self.tabs, padding=12)
        interaction = ttk.Frame(self.tabs, padding=12)
        spatial_graph = ttk.Frame(self.tabs, padding=12)
        dashboard = ttk.Frame(self.tabs, padding=10)
        self.dashboard_tab = dashboard
        self.tabs.add(composition, text="Gene Composition")
        self.tabs.add(enrichment, text="Interface Enrichment")
        self.tabs.add(interaction, text="Cx/Sx Interaction")
        self.tabs.add(spatial_graph, text="Spatial Graph & Neighborhood — Experimental")
        self.tabs.add(dashboard, text="Results Dashboard")
        self._build_composition(composition)
        self._build_enrichment(enrichment)
        self._build_interaction(interaction)
        self._build_spatial_graph(spatial_graph)
        self._build_dashboard(dashboard)

        footer = ttk.Frame(self)
        footer.pack(fill="x", pady=(9, 0))
        self.progress = ttk.Progressbar(footer, mode="indeterminate", length=170)
        self.progress.pack(side="left")
        self.status = ttk.Label(footer, text="Ready. Select one or more h5ad samples on the left.")
        self.status.pack(side="left", padx=10)
        self.open_button = ttk.Button(footer, text="Open latest results", command=self._open_latest, state="disabled")
        self.open_button.pack(side="right")

    @staticmethod
    def _description(parent, title: str, body: str, outputs: str) -> None:
        ttk.Label(parent, text=title, font=("Segoe UI Semibold", 13)).pack(anchor="w")
        ttk.Label(parent, text=body, wraplength=890, justify="left").pack(anchor="w", fill="x", pady=(4, 10))
        box = ttk.LabelFrame(parent, text="Reproducible outputs", padding=9)
        box.pack(fill="x")
        ttk.Label(box, text=outputs, wraplength=870, justify="left").pack(anchor="w", fill="x")

    def _build_composition(self, parent) -> None:
        self._description(
            parent,
            "Within-program gene contribution",
            "Calculates the relative contribution of every requested Cx and Sx gene on the same transformed expression scale used by the core workflow. Missing genes remain visible and are marked explicitly.",
            "gene_composition.csv, 300-dpi PNG, vector PDF, analysis_metadata.json, and a run manifest.",
        )
        ttk.Button(parent, text="Run Gene Composition", style="Primary.TButton", command=lambda: self._run("composition")).pack(fill="x", pady=(14, 0))

    def _build_enrichment(self, parent) -> None:
        self._description(
            parent,
            "Interface-like versus non-interface composition",
            "Uses the unchanged core interface call. Reports group means, composition percentages, fold enrichment, Hedges' g, two-sided Mann-Whitney p-values, and Benjamini-Hochberg FDR.",
            "interface_enrichment.csv, 300-dpi PNG, vector PDF, analysis_metadata.json, and a run manifest.",
        )
        ttk.Label(
            parent,
            text="If either region has too few observations, descriptive values are retained and inferential statistics are reported as unavailable.",
            foreground="#7c2d12",
            wraplength=870,
        ).pack(anchor="w", pady=(10, 0))
        ttk.Button(parent, text="Run Interface Enrichment", style="Primary.TButton", command=lambda: self._run("enrichment")).pack(fill="x", pady=(14, 0))

    def _build_interaction(self, parent) -> None:
        self._description(
            parent,
            "Local Cx/Sx spatial interaction",
            "Builds the same six-neighbor spatial graph and quantifies neighborhood coexistence, antagonism, balance, weighted overlap, and Cx/Sx boundary mixing. No correlation is used as the interaction result.",
            "interaction_summary.csv, per-spot metrics CSV, 300-dpi PNG, vector PDF, analysis_metadata.json, and a run manifest.",
        )
        settings = ttk.LabelFrame(parent, text="Spatial null model", padding=9)
        settings.pack(fill="x", pady=(11, 0))
        ttk.Label(settings, text="Permutations").grid(row=0, column=0, sticky="w")
        ttk.Spinbox(settings, from_=1, to=9999, textvariable=self.permutations_var, width=10).grid(row=0, column=1, sticky="w", padx=(6, 18))
        ttk.Label(settings, text="Seed").grid(row=0, column=2, sticky="w")
        ttk.Spinbox(settings, from_=0, to=2147483647, textvariable=self.seed_var, width=14).grid(row=0, column=3, sticky="w", padx=6)
        ttk.Label(
            settings,
            text="Sx values are permuted across fixed coordinates, then local neighborhoods are recomputed.",
            foreground="#4b5563",
        ).grid(row=1, column=0, columnspan=4, sticky="w", pady=(7, 0))
        ttk.Button(parent, text="Run Cx/Sx Interaction", style="Primary.TButton", command=lambda: self._run("interaction")).pack(fill="x", pady=(14, 0))

    def _build_spatial_graph(self, parent) -> None:
        canvas = tk.Canvas(parent, highlightthickness=0)
        scrollbar = ttk.Scrollbar(parent, orient="vertical", command=canvas.yview)
        content = ttk.Frame(canvas)
        window = canvas.create_window((0, 0), window=content, anchor="nw")
        content.bind("<Configure>", lambda _event: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.bind("<Configure>", lambda event: canvas.itemconfigure(window, width=event.width))
        canvas.bind("<MouseWheel>", lambda event: canvas.yview_scroll(int(-event.delta / 120), "units"))
        canvas.configure(yscrollcommand=scrollbar.set)
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        parent = content
        ttk.Label(parent, text="Spatial Graph & Neighborhood — Experimental", font=("Segoe UI Semibold", 13)).pack(anchor="w")
        ttk.Label(
            parent,
            text=(
                "These analyses provide exploratory spatial association and organization summaries. "
                "They do not establish causal, physical, or biological cell-cell interactions."
            ),
            foreground="#b45309",
            wraplength=890,
            justify="left",
        ).pack(anchor="w", fill="x", pady=(4, 6))
        ttk.Label(
            parent,
            text=(
                "v0.4 optional context module. It builds a sparse spatial graph, computes graph QC, optional H_expr and V_expr "
                "context fields, neighborhood enrichment, binary-mask association, and continuous edge statistics. "
                "It does not alter R(x), Type A/B/C calls, or existing Main Mapper outputs."
            ),
            wraplength=890,
            justify="left",
        ).pack(anchor="w", fill="x", pady=(4, 10))

        audit_box = ttk.LabelFrame(parent, text="1  Input audit", padding=9)
        audit_box.pack(fill="x")
        ttk.Label(
            audit_box,
            text="Every sample is audited without altering AnnData.X. JSON/CSV record preprocessing guess, platform, matrix storage, coordinates, tissue metadata, duplicate names, library-size summaries, and warnings.",
            foreground="#4b5563", wraplength=850, justify="left",
        ).pack(anchor="w")

        graph_box = ttk.LabelFrame(parent, text="2  Graph setup", padding=9)
        graph_box.pack(fill="x")
        row = ttk.Frame(graph_box); row.pack(fill="x")
        ttk.Label(row, text="Method").grid(row=0, column=0, sticky="w")
        ttk.Combobox(row, textvariable=self.graph_method_var, values=("radius", "lattice", "knn"), state="readonly", width=12).grid(row=1, column=0, sticky="ew")
        ttk.Label(row, text="Radius").grid(row=0, column=1, sticky="w", padx=(8, 0))
        ttk.Entry(row, textvariable=self.graph_radius_var, width=10).grid(row=1, column=1, sticky="ew", padx=(8, 0))
        ttk.Label(row, text="K").grid(row=0, column=2, sticky="w", padx=(8, 0))
        ttk.Spinbox(row, from_=1, to=50, textvariable=self.graph_k_var, width=7).grid(row=1, column=2, sticky="ew", padx=(8, 0))
        ttk.Label(row, text="Weighting").grid(row=0, column=3, sticky="w", padx=(8, 0))
        ttk.Combobox(row, textvariable=self.graph_weighting_var, values=("binary", "inverse_distance", "gaussian"), state="readonly", width=16).grid(row=1, column=3, sticky="ew", padx=(8, 0))
        ttk.Label(row, text="KNN sym.").grid(row=0, column=4, sticky="w", padx=(8, 0))
        ttk.Combobox(row, textvariable=self.graph_sym_var, values=("union", "mutual"), state="readonly", width=10).grid(row=1, column=4, sticky="ew", padx=(8, 0))
        for column in range(5):
            row.columnconfigure(column, weight=1)
        calibration = ttk.Frame(graph_box); calibration.pack(fill="x", pady=(7, 0))
        ttk.Label(calibration, text="Coordinate unit").grid(row=0, column=0, sticky="w")
        ttk.Combobox(calibration, textvariable=self.graph_coordinate_unit_var, values=("native", "pixel", "micrometer"), state="readonly", width=13).grid(row=1, column=0, sticky="ew")
        ttk.Label(calibration, text="µm per coordinate unit").grid(row=0, column=1, sticky="w", padx=(8, 0))
        ttk.Entry(calibration, textvariable=self.graph_coordinate_scale_var, width=12).grid(row=1, column=1, sticky="ew", padx=(8, 0))
        ttk.Label(calibration, text="Scale source").grid(row=0, column=2, sticky="w", padx=(8, 0))
        ttk.Entry(calibration, textvariable=self.graph_scale_source_var).grid(row=1, column=2, sticky="ew", padx=(8, 0))
        calibration.columnconfigure(2, weight=1)
        self.graph_coordinate_status = ttk.Label(graph_box, text="Coordinate unit: native", foreground="#4b5563")
        self.graph_coordinate_status.pack(anchor="w", pady=(5, 0))
        ttk.Label(
            graph_box,
            text="Radius graph is the default. It is a calibrated physical-radius graph only when a valid micrometer scale and scale source are recorded. KNN remains an auxiliary robustness graph.",
            foreground="#4b5563",
            wraplength=850,
        ).pack(anchor="w", pady=(7, 0))

        graph_qc_box = ttk.LabelFrame(parent, text="3  Graph QC", padding=9)
        graph_qc_box.pack(fill="x", pady=(9, 0))
        ttk.Label(graph_qc_box, text="Exports degree, isolated fraction, connected components, edge-distance quantiles, duplicate coordinates, graph density, long-edge warnings, and requested/effective graph provenance.", foreground="#4b5563", wraplength=850).pack(anchor="w")

        context = ttk.LabelFrame(parent, text="4  Context fields", padding=9)
        context.pack(fill="x", pady=(9, 0))
        top = ttk.Frame(context); top.pack(fill="x")
        ttk.Checkbutton(top, text="Enable H_expr hypoxia-associated expression field", variable=self.graph_enable_h_var).grid(row=0, column=0, sticky="w")
        ttk.Checkbutton(top, text="Enable V_expr endothelial/angiogenic expression proxy", variable=self.graph_enable_v_var).grid(row=0, column=1, sticky="w", padx=(10, 0))
        ttk.Label(top, text="Smoothing").grid(row=0, column=2, sticky="w", padx=(10, 0))
        ttk.Combobox(top, textvariable=self.graph_smoothing_var, values=("none", "graph_mean"), state="readonly", width=12).grid(row=0, column=3, sticky="w", padx=(5, 0))
        ttk.Label(
            context,
            text=(
                "H_expr is a hypoxia-associated expression field. V_expr is an endothelial/angiogenic expression proxy. "
                "V_expr is not a direct measure of vessel density, perfusion, or functional blood supply."
            ),
            foreground="#4b5563",
            wraplength=850,
            justify="left",
        ).pack(anchor="w", pady=(7, 0))
        ttk.Label(
            context,
            text=(
                "Graph-smoothed context fields are intended for visualization and exploratory sensitivity analysis. "
                "Association statistics computed on fields smoothed over the same graph may be inflated and should not "
                "be interpreted as independent confirmatory evidence."
            ),
            foreground="#b45309",
            wraplength=850,
            justify="left",
        ).pack(anchor="w", pady=(5, 0))
        methods = ttk.Frame(context); methods.pack(fill="x", pady=(8, 0))
        ttk.Label(methods, text="H method").grid(row=0, column=0, sticky="w")
        ttk.Combobox(methods, textvariable=self.graph_h_method_var, values=("raw_mean", "z_score_mean", "rank_quantile"), state="readonly", width=15).grid(row=1, column=0, sticky="ew")
        ttk.Label(methods, text="H-high q").grid(row=0, column=1, sticky="w", padx=(8, 0))
        ttk.Entry(methods, textvariable=self.graph_h_quantile_var, width=8).grid(row=1, column=1, sticky="ew", padx=(8, 0))
        ttk.Label(methods, text="V method").grid(row=0, column=2, sticky="w", padx=(8, 0))
        ttk.Combobox(methods, textvariable=self.graph_v_method_var, values=("raw_mean", "z_score_mean", "rank_quantile"), state="readonly", width=15).grid(row=1, column=2, sticky="ew", padx=(8, 0))
        ttk.Label(methods, text="V-high q").grid(row=0, column=3, sticky="w", padx=(8, 0))
        ttk.Entry(methods, textvariable=self.graph_v_quantile_var, width=8).grid(row=1, column=3, sticky="ew", padx=(8, 0))
        ttk.Label(methods, text="Min coverage").grid(row=0, column=4, sticky="w", padx=(8, 0))
        ttk.Entry(methods, textvariable=self.graph_min_coverage_var, width=8).grid(row=1, column=4, sticky="ew", padx=(8, 0))
        ttk.Label(methods, text="Min spot fraction").grid(row=0, column=5, sticky="w", padx=(8, 0))
        ttk.Entry(methods, textvariable=self.graph_min_spot_fraction_var, width=8).grid(row=1, column=5, sticky="ew", padx=(8, 0))
        ttk.Checkbutton(methods, text="Allow low coverage", variable=self.graph_allow_low_coverage_var).grid(row=1, column=6, sticky="w", padx=(8, 0))
        for column in range(7):
            methods.columnconfigure(column, weight=1)
        gene_rows = ttk.Frame(context); gene_rows.pack(fill="x", pady=(8, 0))
        ttk.Label(gene_rows, text="H genes").grid(row=0, column=0, sticky="nw")
        ttk.Entry(gene_rows, textvariable=self.graph_h_genes_var).grid(row=0, column=1, sticky="ew", padx=(7, 0))
        ttk.Label(gene_rows, text="V genes").grid(row=1, column=0, sticky="nw", pady=(5, 0))
        ttk.Entry(gene_rows, textvariable=self.graph_v_genes_var).grid(row=1, column=1, sticky="ew", padx=(7, 0), pady=(5, 0))
        gene_rows.columnconfigure(1, weight=1)

        overlap = ttk.LabelFrame(parent, text="5  Same-spot overlap", padding=9); overlap.pack(fill="x", pady=(9, 0))
        ttk.Label(overlap, text="Optional binary mask A column").grid(row=0, column=0, sticky="w")
        ttk.Entry(overlap, textvariable=self.graph_mask_a_var).grid(row=1, column=0, sticky="ew")
        ttk.Label(overlap, text="Optional binary mask B column").grid(row=0, column=1, sticky="w", padx=(8, 0))
        ttk.Entry(overlap, textvariable=self.graph_mask_b_var).grid(row=1, column=1, sticky="ew", padx=(8, 0))
        overlap.columnconfigure(0, weight=1); overlap.columnconfigure(1, weight=1)

        neighborhood = ttk.LabelFrame(parent, text="6  Neighborhood enrichment", padding=9); neighborhood.pack(fill="x", pady=(9, 0))
        ttk.Label(neighborhood, text="Optional adata.obs categorical state column (blank = SpatialTX spot-level states)").grid(row=0, column=0, sticky="w")
        ttk.Entry(neighborhood, textvariable=self.graph_label_source_var).grid(row=1, column=0, sticky="ew")
        ttk.Combobox(neighborhood, textvariable=self.graph_label_mode_var, values=("categorical_state", "binary_mask"), state="readonly", width=18).grid(row=1, column=1, sticky="ew", padx=(8, 0))
        neighborhood.columnconfigure(0, weight=1)

        continuous = ttk.LabelFrame(parent, text="7  Continuous edge association", padding=9); continuous.pack(fill="x", pady=(9, 0))
        ttk.Label(continuous, text="Optional X column").grid(row=0, column=0, sticky="w")
        ttk.Entry(continuous, textvariable=self.graph_x_column_var).grid(row=1, column=0, sticky="ew")
        ttk.Combobox(continuous, textvariable=self.graph_x_mode_var, values=("continuous_score", "proportion_composition"), state="readonly", width=21).grid(row=1, column=1, padx=(6, 0))
        ttk.Label(continuous, text="Optional Y column").grid(row=0, column=2, sticky="w", padx=(8, 0))
        ttk.Entry(continuous, textvariable=self.graph_y_column_var).grid(row=1, column=2, sticky="ew", padx=(8, 0))
        ttk.Combobox(continuous, textvariable=self.graph_y_mode_var, values=("continuous_score", "proportion_composition"), state="readonly", width=21).grid(row=1, column=3, padx=(6, 0))
        continuous.columnconfigure(0, weight=1); continuous.columnconfigure(2, weight=1)

        analysis = ttk.LabelFrame(parent, text="8  Permutation scope", padding=9); analysis.pack(fill="x", pady=(9, 0))
        ttk.Label(analysis, text="Scope").grid(row=0, column=0, sticky="w")
        ttk.Combobox(analysis, textvariable=self.graph_permutation_scope_var, values=("whole_slide", "within_connected_components", "within_user_strata"), state="readonly", width=28).grid(row=1, column=0, sticky="ew")
        ttk.Label(analysis, text="Stratification column").grid(row=0, column=1, sticky="w", padx=(8, 0))
        ttk.Entry(analysis, textvariable=self.graph_strata_column_var).grid(row=1, column=1, sticky="ew", padx=(8, 0))
        ttk.Label(analysis, text="Permutations").grid(row=0, column=2, sticky="w", padx=(8, 0))
        ttk.Spinbox(analysis, from_=1, to=99999, textvariable=self.graph_permutations_var, width=10).grid(row=1, column=2, sticky="ew", padx=(8, 0))
        ttk.Label(analysis, text="Seed").grid(row=0, column=3, sticky="w", padx=(8, 0))
        ttk.Spinbox(analysis, from_=0, to=2147483647, textvariable=self.graph_seed_var, width=14).grid(row=1, column=3, sticky="ew", padx=(8, 0))
        ttk.Checkbutton(analysis, text="Restrict to in_tissue=1 when available", variable=self.graph_tissue_only_var).grid(row=2, column=0, columnspan=2, sticky="w", pady=(6, 0))
        ttk.Label(analysis, text="Disconnected tissue fragments: within_connected_components is recommended. Permutations never cross samples.", foreground="#4b5563").grid(row=2, column=2, columnspan=2, sticky="w", pady=(6, 0))
        analysis.columnconfigure(0, weight=1); analysis.columnconfigure(1, weight=1)

        robustness = ttk.LabelFrame(parent, text="9  Robustness comparison", padding=9); robustness.pack(fill="x", pady=(9, 0))
        ttk.Checkbutton(robustness, text="Compare radius, lattice, and KNN association direction (optional)", variable=self.graph_robustness_var).pack(anchor="w")
        ttk.Label(robustness, text="The comparison is supplementary and never replaces the primary graph result.", foreground="#4b5563").pack(anchor="w", pady=(4, 0))

        actions = ttk.Frame(parent)
        actions.pack(fill="x", pady=(10, 0))
        ttk.Button(actions, text="Run Spatial Graph & Neighborhood", style="Primary.TButton", command=self._run_spatial_graph).pack(side="left", fill="x", expand=True)
        self.graph_cancel_button = ttk.Button(actions, text="Cancel analysis", command=self._cancel_spatial_graph, state="disabled")
        self.graph_cancel_button.pack(side="left", padx=(8, 0))

        preview = ttk.LabelFrame(parent, text="10  Export preview", padding=8)
        preview.pack(fill="both", expand=True, pady=(9, 0))
        self.graph_preview = tk.Text(preview, height=8, wrap="word", state="disabled", background="#f8fafc")
        self.graph_preview.pack(fill="both", expand=True)

    def _build_dashboard(self, parent) -> None:
        heading = ttk.Frame(parent)
        heading.pack(fill="x")
        ttk.Label(heading, text="Advanced Analysis Results", font=("Segoe UI Semibold", 13)).pack(side="left")
        ttk.Button(heading, text="Clear dashboard", command=self._clear_dashboard).pack(side="right")
        ttk.Label(
            parent,
            text="The upper table summarizes every completed module and sample. Select a row to inspect its complete CSV table below.",
            foreground="#4b5563",
        ).pack(anchor="w", pady=(3, 8))

        overview_box = ttk.LabelFrame(parent, text="At-a-glance summary", padding=7)
        overview_box.pack(fill="x")
        columns = ("module", "sample", "status", "metric_1", "metric_2", "metric_3")
        self.overview_tree = ttk.Treeview(overview_box, columns=columns, show="headings", height=7, selectmode="browse")
        headings = {
            "module": "Module", "sample": "Sample", "status": "Status",
            "metric_1": "Key result 1", "metric_2": "Key result 2", "metric_3": "Key result 3",
        }
        widths = {"module": 150, "sample": 150, "status": 70, "metric_1": 220, "metric_2": 220, "metric_3": 220}
        for column in columns:
            self.overview_tree.heading(column, text=headings[column])
            self.overview_tree.column(column, width=widths[column], minwidth=70, anchor="w")
        overview_x = ttk.Scrollbar(overview_box, orient="horizontal", command=self.overview_tree.xview)
        overview_y = ttk.Scrollbar(overview_box, orient="vertical", command=self.overview_tree.yview)
        self.overview_tree.configure(xscrollcommand=overview_x.set, yscrollcommand=overview_y.set)
        self.overview_tree.grid(row=0, column=0, sticky="nsew")
        overview_y.grid(row=0, column=1, sticky="ns")
        overview_x.grid(row=1, column=0, sticky="ew")
        overview_box.columnconfigure(0, weight=1)
        overview_box.rowconfigure(0, weight=1)
        self.overview_tree.tag_configure("error", foreground="#b91c1c")
        self.overview_tree.bind("<<TreeviewSelect>>", self._show_dashboard_detail)

        detail_header = ttk.Frame(parent)
        detail_header.pack(fill="x", pady=(9, 4))
        self.detail_label = ttk.Label(detail_header, text="Detailed table — select a summary row", font=("Segoe UI Semibold", 10))
        self.detail_label.pack(side="left")
        self.open_figure_button = ttk.Button(detail_header, text="Open figure", command=self._open_selected_figure, state="disabled")
        self.open_figure_button.pack(side="right")
        self.open_csv_button = ttk.Button(detail_header, text="Open CSV", command=self._open_selected_csv, state="disabled")
        self.open_csv_button.pack(side="right", padx=(0, 5))

        detail_box = ttk.Frame(parent)
        detail_box.pack(fill="both", expand=True)
        self.detail_tree = ttk.Treeview(detail_box, show="headings", height=12)
        detail_x = ttk.Scrollbar(detail_box, orient="horizontal", command=self.detail_tree.xview)
        detail_y = ttk.Scrollbar(detail_box, orient="vertical", command=self.detail_tree.yview)
        self.detail_tree.configure(xscrollcommand=detail_x.set, yscrollcommand=detail_y.set)
        self.detail_tree.grid(row=0, column=0, sticky="nsew")
        detail_y.grid(row=0, column=1, sticky="ns")
        detail_x.grid(row=1, column=0, sticky="ew")
        detail_box.columnconfigure(0, weight=1)
        detail_box.rowconfigure(0, weight=1)
        self.dashboard_note = ttk.Label(parent, text="No Advanced Analysis results loaded yet.", foreground="#4b5563")
        self.dashboard_note.pack(anchor="w", pady=(5, 0))

    @staticmethod
    def _display_value(value) -> str:
        if pd.isna(value):
            return ""
        if isinstance(value, (float, np.floating)):
            magnitude = abs(float(value))
            if magnitude and (magnitude < .001 or magnitude >= 10000):
                return f"{float(value):.3e}"
            return f"{float(value):.4f}".rstrip("0").rstrip(".")
        return str(value)

    @staticmethod
    def _numeric(table: pd.DataFrame, column: str) -> pd.Series:
        if column not in table:
            return pd.Series(np.nan, index=table.index, dtype=float)
        return pd.to_numeric(table[column], errors="coerce")

    def _summary_metrics(self, module: str, table: pd.DataFrame) -> tuple[str, str, str]:
        if table.empty:
            return "No table rows", "", ""
        if module == "composition":
            present = table[table.get("status", "").eq("present")].copy()
            contribution = self._numeric(present, "relative_contribution_percent")
            summaries: list[str] = []
            for program in ("Cx", "Sx"):
                group = present[present.get("program", "").eq(program)]
                values = self._numeric(group, "relative_contribution_percent")
                if len(group) and values.notna().any():
                    index = values.idxmax()
                    summaries.append(f"{program} top: {group.loc[index, 'gene']} ({values.loc[index]:.1f}%)")
                else:
                    summaries.append(f"{program} top: unavailable")
            missing = int(table.get("status", pd.Series(index=table.index, dtype=str)).eq("missing").sum())
            return summaries[0], summaries[1], f"Present {len(present)} / Missing {missing}"
        if module == "enrichment":
            present = table[table.get("status", "").eq("present")].copy()
            n_interface = self._numeric(present, "n_interface").dropna()
            n_noninterface = self._numeric(present, "n_noninterface").dropna()
            significant = int(present.get("significant_fdr_0_05", pd.Series(False, index=present.index)).astype(str).str.lower().eq("true").sum())
            effects = self._numeric(present, "hedges_g").abs()
            strongest = "Strongest effect: unavailable"
            if effects.notna().any():
                index = effects.idxmax()
                strongest = f"Strongest: {present.loc[index, 'gene']} (|g|={effects.loc[index]:.2f})"
            group_sizes = f"Interface {int(n_interface.iloc[0]) if len(n_interface) else 0} / Other {int(n_noninterface.iloc[0]) if len(n_noninterface) else 0}"
            return group_sizes, f"FDR < 0.05: {significant}/{len(present)} genes", strongest
        row = table.iloc[0]
        return (
            f"Coexistence {self._display_value(row.get('coexistence_index'))}",
            f"Antagonism {self._display_value(row.get('antagonism_index'))}",
            f"Overlap {self._display_value(row.get('spatial_overlap_index'))} / Balance {self._display_value(row.get('balance_index'))}",
        )

    def _add_dashboard_run(self, module: str, run_dir: Path, manifest: pd.DataFrame) -> None:
        for _, row in manifest.iterrows():
            sample = str(row.get("sample", ""))
            for item_id, record in list(self.dashboard_records.items()):
                if record["module"] == module and record["sample"] == sample:
                    self.overview_tree.delete(item_id)
                    del self.dashboard_records[item_id]
            status = str(row.get("status", ""))
            csv_path = Path(str(row.get("table_csv", ""))) if status == "ok" else None
            figure_path = Path(str(row.get("figure_png", ""))) if status == "ok" else None
            table = pd.DataFrame()
            metrics = ("", "", "")
            if csv_path and csv_path.is_file():
                try:
                    table = pd.read_csv(csv_path)
                    metrics = self._summary_metrics(module, table)
                except Exception as exc:
                    status = f"error: cannot read table ({exc})"
            self.dashboard_counter += 1
            item_id = f"advanced_result_{self.dashboard_counter}"
            self.overview_tree.insert(
                "", "end", iid=item_id,
                values=(MODULE_LABELS[module], sample, "OK" if status == "ok" else "ERROR", *metrics),
                tags=() if status == "ok" else ("error",),
            )
            self.dashboard_records[item_id] = {
                "module": module, "sample": sample, "status": status, "table": table,
                "csv": csv_path, "figure": figure_path, "run_dir": run_dir,
            }
        self.dashboard_note.configure(text=f"Loaded {len(self.dashboard_records)} module/sample result(s). Detailed tables show up to 2,000 rows.")
        children = self.overview_tree.get_children()
        if children:
            self.overview_tree.selection_set(children[0])
            self.overview_tree.focus(children[0])
            self._show_dashboard_detail()

    def _show_dashboard_detail(self, _event=None) -> None:
        selection = self.overview_tree.selection()
        if not selection:
            return
        record = self.dashboard_records.get(selection[0])
        if not record:
            return
        self.selected_record = record
        table: pd.DataFrame = record["table"]
        self.detail_tree.delete(*self.detail_tree.get_children())
        columns = [str(column) for column in table.columns]
        self.detail_tree.configure(columns=columns)
        for column in columns:
            values = [self._display_value(value) for value in table[column].head(30)]
            width = min(260, max(90, 8 * max([len(column), *(len(value) for value in values)])))
            self.detail_tree.heading(column, text=column)
            self.detail_tree.column(column, width=width, minwidth=70, anchor="w", stretch=False)
        for index, row in table.head(2000).iterrows():
            self.detail_tree.insert("", "end", iid=f"detail_{index}", values=[self._display_value(row[column]) for column in table.columns])
        self.detail_label.configure(text=f"Detailed table — {MODULE_LABELS[record['module']]} / {record['sample']} ({len(table)} rows)")
        self.open_csv_button.configure(state="normal" if record["csv"] and record["csv"].is_file() else "disabled")
        self.open_figure_button.configure(state="normal" if record["figure"] and record["figure"].is_file() else "disabled")

    def _clear_dashboard(self) -> None:
        self.overview_tree.delete(*self.overview_tree.get_children())
        self.detail_tree.delete(*self.detail_tree.get_children())
        self.detail_tree.configure(columns=())
        self.dashboard_records.clear()
        self.selected_record = None
        self.detail_label.configure(text="Detailed table — select a summary row")
        self.dashboard_note.configure(text="No Advanced Analysis results loaded yet.")
        self.open_csv_button.configure(state="disabled")
        self.open_figure_button.configure(state="disabled")

    @staticmethod
    def _open_path(path: Path) -> None:
        if os.name == "nt":
            os.startfile(path)  # type: ignore[attr-defined]
        else:
            subprocess.Popen(["open" if os.uname().sysname == "Darwin" else "xdg-open", str(path)])

    def _open_selected_csv(self) -> None:
        if self.selected_record and self.selected_record["csv"] and self.selected_record["csv"].is_file():
            try:
                self._open_path(self.selected_record["csv"])
            except Exception as exc:
                messagebox.showerror("Open CSV", str(exc), parent=self)

    def _open_selected_figure(self) -> None:
        if self.selected_record and self.selected_record["figure"] and self.selected_record["figure"].is_file():
            try:
                self._open_path(self.selected_record["figure"])
            except Exception as exc:
                messagebox.showerror("Open figure", str(exc), parent=self)

    def _browse_output(self) -> None:
        value = filedialog.askdirectory(title="Advanced Analysis output root", parent=self)
        if value:
            self.output_var.set(value)

    def _append_graph_preview(self, text: str) -> None:
        self.graph_preview.configure(state="normal")
        self.graph_preview.insert("end", text.rstrip() + "\n")
        self.graph_preview.see("end")
        self.graph_preview.configure(state="disabled")

    def _update_coordinate_status(self, *_args) -> None:
        if not hasattr(self, "graph_coordinate_status"):
            return
        unit = self.graph_coordinate_unit_var.get().strip() or "native"
        scale = self.graph_coordinate_scale_var.get().strip()
        source = self.graph_scale_source_var.get().strip()
        if unit == "micrometer" and source:
            text = "Coordinate unit: micrometer, calibrated"
        elif unit in {"native", "pixel"} and scale and source:
            text = f"Coordinate unit: {unit}; calibrated to micrometer"
        elif unit == "micrometer":
            text = "Coordinate unit: micrometer, calibration metadata incomplete"
        else:
            text = f"Coordinate unit: {unit}"
        self.graph_coordinate_status.configure(text=text)

    def _build_graph_config(self):
        from .graph.builder import GraphBuildConfig
        from .graph.runner import SpatialGraphAnalysisConfig

        radius_text = self.graph_radius_var.get().strip()
        radius = float(radius_text) if radius_text else None
        scale_text = self.graph_coordinate_scale_var.get().strip()
        coordinate_scale = float(scale_text) if scale_text else None
        label_source = self.graph_label_source_var.get().strip() or "auto_spatialtx_states"
        return SpatialGraphAnalysisConfig(
            graph=GraphBuildConfig(
                method=self.graph_method_var.get(),
                radius=radius,
                k=int(self.graph_k_var.get()),
                weighting=self.graph_weighting_var.get(),
                symmetrization=self.graph_sym_var.get(),
                random_seed=int(self.graph_seed_var.get()),
                coordinate_unit=self.graph_coordinate_unit_var.get(),
                coordinate_scale=coordinate_scale,
                scale_source=self.graph_scale_source_var.get().strip(),
            ),
            enable_h=bool(self.graph_enable_h_var.get()),
            enable_v=bool(self.graph_enable_v_var.get()),
            h_genes=parse_gene_text(self.graph_h_genes_var.get()),
            v_genes=parse_gene_text(self.graph_v_genes_var.get()),
            h_score_method=self.graph_h_method_var.get(),
            v_score_method=self.graph_v_method_var.get(),
            h_high_quantile=float(self.graph_h_quantile_var.get()),
            v_high_quantile=float(self.graph_v_quantile_var.get()),
            context_smoothing=self.graph_smoothing_var.get(),
            min_gene_coverage=float(self.graph_min_coverage_var.get()),
            context_min_spot_fraction=float(self.graph_min_spot_fraction_var.get()),
            allow_low_coverage=bool(self.graph_allow_low_coverage_var.get()),
            label_source=label_source,
            label_mode=self.graph_label_mode_var.get(),
            user_mask_a_column=self.graph_mask_a_var.get().strip(),
            user_mask_b_column=self.graph_mask_b_var.get().strip(),
            continuous_x_column=self.graph_x_column_var.get().strip(),
            continuous_y_column=self.graph_y_column_var.get().strip(),
            continuous_x_mode=self.graph_x_mode_var.get(),
            continuous_y_mode=self.graph_y_mode_var.get(),
            permutation_scope=self.graph_permutation_scope_var.get(),
            stratification_column=self.graph_strata_column_var.get().strip(),
            tissue_only_restriction=bool(self.graph_tissue_only_var.get()),
            run_graph_robustness=bool(self.graph_robustness_var.get()),
            permutations=int(self.graph_permutations_var.get()),
            seed=int(self.graph_seed_var.get()),
        )

    def _set_enabled(self, enabled: bool) -> None:
        for child in self.winfo_children():
            self._set_widget_enabled(child, enabled)

    def _set_widget_enabled(self, widget, enabled: bool) -> None:
        if isinstance(widget, (ttk.Button, ttk.Entry, ttk.Spinbox, ttk.Combobox, ttk.Checkbutton, ttk.Notebook)):
            widget.state(["!disabled"] if enabled else ["disabled"])
        for child in widget.winfo_children():
            self._set_widget_enabled(child, enabled)

    def _run_spatial_graph(self) -> None:
        if self.busy:
            return
        try:
            samples = self.get_samples()
            if not samples:
                raise ValueError("Select at least one scanned h5ad sample in the main sample list.")
            c_genes, s_genes = self.get_genes()
            c_q, s_q, g_q = self.get_quantiles()
            output = self.output_var.get().strip()
            if not output:
                raise ValueError("Choose an output root.")
            config = self._build_graph_config()
            config.c_q, config.s_q, config.g_q = c_q, s_q, g_q
        except Exception as exc:
            messagebox.showwarning("Spatial Graph & Neighborhood", str(exc), parent=self)
            return
        self.busy = True
        self.graph_cancel_event = Event()
        self._set_enabled(False)
        self.graph_cancel_button.state(["!disabled"])
        self.progress.start(12)
        self.status.configure(text="Running Spatial Graph & Neighborhood...")
        self.graph_preview.configure(state="normal")
        self.graph_preview.delete("1.0", "end")
        self.graph_preview.configure(state="disabled")
        self._append_graph_preview("Started Spatial Graph and Neighborhood Analysis.")

        def progress(message: str) -> None:
            self.events.put(("graph_progress", message))

        def work() -> None:
            try:
                from .graph.runner import run_spatial_graph_neighborhood_batch

                run_dir, manifest = run_spatial_graph_neighborhood_batch(
                    samples,
                    output,
                    c_genes,
                    s_genes,
                    config,
                    progress,
                    self.graph_cancel_event,
                )
                self.events.put(("graph_done", run_dir, manifest))
            except Exception as exc:
                self.events.put(("graph_error", exc))

        threading.Thread(target=work, daemon=True).start()

    def _cancel_spatial_graph(self) -> None:
        if self.graph_cancel_event is not None:
            self.graph_cancel_event.set()
            self.status.configure(text="Cancellation requested; finishing current sample if already running...")
            self._append_graph_preview("Cancellation requested.")

    def _run(self, module: str) -> None:
        if self.busy:
            return
        try:
            samples = self.get_samples()
            if not samples:
                raise ValueError("Select at least one scanned h5ad sample in the main sample list.")
            c_genes, s_genes = self.get_genes()
            c_q, s_q, g_q = self.get_quantiles()
            output = self.output_var.get().strip()
            if not output:
                raise ValueError("Choose an output root.")
            permutations = int(self.permutations_var.get())
            seed = int(self.seed_var.get())
        except Exception as exc:
            messagebox.showwarning("Advanced Analysis", str(exc), parent=self)
            return
        self.busy = True
        self._set_enabled(False)
        self.progress.start(12)
        modules = ["composition", "enrichment", "interaction"] if module == "all" else [module]
        if module == "all":
            self._clear_dashboard()
        label = "all three Advanced Analysis modules" if module == "all" else MODULE_LABELS[module]
        self.status.configure(text=f"Running {label}...")

        def progress(message: str) -> None:
            self.events.put(("progress", message))

        def work() -> None:
            try:
                results = []
                for current in modules:
                    progress(f"Starting {MODULE_LABELS[current]}...")
                    result = run_advanced_batch(
                        current, samples, output, c_genes, s_genes, progress, c_q, s_q, g_q, permutations, seed
                    )
                    results.append((current, result))
                self.events.put(("done", module, results))
            except Exception as exc:
                self.events.put(("error", module, exc))

        threading.Thread(target=work, daemon=True).start()

    def _poll(self) -> None:
        try:
            while True:
                event = self.events.get_nowait()
                if event[0] == "progress":
                    self.status.configure(text=event[1])
                elif event[0] == "graph_progress":
                    self.status.configure(text=event[1])
                    self._append_graph_preview(event[1])
                elif event[0] == "done":
                    _, module, results = event
                    self.busy = False; self._set_enabled(True); self.progress.stop()
                    successes = total = 0
                    for current, (run_dir, manifest) in results:
                        self.last_run = run_dir
                        self._add_dashboard_run(current, run_dir, manifest)
                        successes += int((manifest["status"] == "ok").sum())
                        total += len(manifest)
                    self.open_button.configure(state="normal")
                    self.tabs.select(self.dashboard_tab)
                    label = "All three modules" if module == "all" else MODULE_LABELS[module]
                    self.status.configure(text=f"{label} complete: {successes}/{total} result(s). Dashboard updated.")
                    if successes != total:
                        messagebox.showwarning("Advanced Analysis", f"Completed with {total - successes} failed result(s). See the dashboard and run manifests.", parent=self)
                elif event[0] == "graph_done":
                    _, run_dir, manifest = event
                    self.busy = False; self._set_enabled(True); self.progress.stop()
                    self.graph_cancel_button.state(["disabled"])
                    self.graph_cancel_event = None
                    self.last_run = run_dir
                    self.open_button.configure(state="normal")
                    ok = int((manifest["status"] == "ok").sum()) if "status" in manifest else 0
                    total = len(manifest)
                    self.status.configure(text=f"Spatial Graph & Neighborhood complete: {ok}/{total} sample(s).")
                    self._append_graph_preview(f"Output folder: {run_dir}")
                    self._append_graph_preview(manifest.to_string(index=False))
                    if ok != total:
                        messagebox.showwarning("Spatial Graph & Neighborhood", f"Completed with {total - ok} failed or cancelled sample(s). See the manifest.", parent=self)
                elif event[0] == "graph_error":
                    _, exc = event
                    self.busy = False; self._set_enabled(True); self.progress.stop()
                    self.graph_cancel_button.state(["disabled"])
                    self.graph_cancel_event = None
                    self.status.configure(text="Spatial Graph & Neighborhood failed.")
                    self._append_graph_preview(f"Error: {exc}")
                    messagebox.showerror("Spatial Graph & Neighborhood", str(exc), parent=self)
                else:
                    _, module, exc = event
                    self.busy = False; self._set_enabled(True); self.progress.stop()
                    label = "All three modules" if module == "all" else MODULE_LABELS[module]
                    self.status.configure(text=f"{label} failed.")
                    self.open_csv_button.configure(state="normal" if self.selected_record and self.selected_record.get("csv") else "disabled")
                    self.open_figure_button.configure(state="normal" if self.selected_record and self.selected_record.get("figure") else "disabled")
                    messagebox.showerror("Advanced Analysis", str(exc), parent=self)
        except queue.Empty:
            pass
        if self.winfo_exists():
            self.after(100, self._poll)

    def _open_latest(self) -> None:
        if not self.last_run or not self.last_run.exists():
            return
        try:
            self._open_path(self.last_run)
        except Exception as exc:
            messagebox.showerror("Open results", str(exc), parent=self)
