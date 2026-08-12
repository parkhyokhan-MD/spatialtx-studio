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

from .comparative.models import ComparativeConfig, ComparativeRunResult, SampleRecord
from .comparative.runner import run_comparative_analysis
from .comparative.validation import load_comparative_manifest, preflight_records, validate_record_structure
from .multi_pair_ui import MultiPairAnalysisPanel


def _sample_display_labels(sample_paths: list[Path]) -> tuple[list[str], dict[str, str]]:
    """Return filename-first labels while retaining an unambiguous path lookup."""
    paths = [Path(path) for path in sample_paths]
    name_counts: dict[str, int] = {}
    for path in paths:
        key = path.name.casefold()
        name_counts[key] = name_counts.get(key, 0) + 1

    labels: list[str] = []
    path_by_label: dict[str, str] = {}
    used_labels: set[str] = set()
    for path in paths:
        label = path.name
        if name_counts[path.name.casefold()] > 1:
            parent_parts = path.parent.parts
            depth = 1
            while True:
                hint = os.sep.join(parent_parts[-depth:]) if parent_parts else str(path.parent)
                candidate = f"{path.name}  [{hint}]"
                if candidate.casefold() not in used_labels:
                    label = candidate
                    break
                if depth < len(parent_parts):
                    depth += 1
                    continue
                suffix = 2
                while f"{candidate} #{suffix}".casefold() in used_labels:
                    suffix += 1
                label = f"{candidate} #{suffix}"
                break
        if label.casefold() in used_labels:
            base_label = label
            suffix = 2
            while f"{base_label} #{suffix}".casefold() in used_labels:
                suffix += 1
            label = f"{base_label} #{suffix}"
        labels.append(label)
        path_by_label[label] = str(path)
        used_labels.add(label.casefold())
    return labels, path_by_label


class ComparativeAnalysisPanel(ttk.Frame):
    """Responsive GUI wrapper around the shared comparative runner."""

    def __init__(
        self,
        parent,
        *,
        get_samples: Callable[[], list[Path]],
        get_genes: Callable[[], tuple[list[str], list[str]]],
        get_quantiles: Callable[[], tuple[float, float, float]],
        get_scoring_options: Callable[[], object],
        get_output: Callable[[], str],
    ) -> None:
        super().__init__(parent)
        self.get_samples = get_samples
        self.get_genes = get_genes
        self.get_quantiles = get_quantiles
        self.get_scoring_options = get_scoring_options
        self.events: queue.Queue[tuple] = queue.Queue()
        self.cancel_event: Event | None = None
        self.busy = False
        self.last_result: ComparativeRunResult | None = None
        self.current_photo = None
        self.figure_label_to_path: dict[str, Path] = {}
        self._resize_job: str | None = None

        self.mode_var = tk.StringVar(value="pairwise")
        self.sample_a_var = tk.StringVar()
        self.sample_b_var = tk.StringVar()
        self.sample_a_display_var = tk.StringVar()
        self.sample_b_display_var = tk.StringVar()
        self._sample_path_by_label: dict[str, str] = {}
        self._sample_label_by_path: dict[str, str] = {}
        self.manifest_var = tk.StringVar()
        self.reference_var = tk.StringVar(value="A")
        self.target_var = tk.StringVar(value="B")
        self.output_var = tk.StringVar(value=get_output())
        self.graph_method_var = tk.StringVar(value="knn")
        self.graph_k_var = tk.IntVar(value=6)
        self.graph_radius_var = tk.StringVar()
        self.test_var = tk.StringVar(value="auto")
        self.seed_var = tk.IntVar(value=42)
        self.h_var = tk.BooleanVar(value=True)
        self.v_var = tk.BooleanVar(value=True)
        self.context_smoothing_var = tk.StringVar(value="none")
        self.figure_var = tk.StringVar()
        self.status_var = tk.StringVar(value="Ready. Comparative outputs are exploratory and non-diagnostic.")
        self.scale_banner_var = tk.StringVar(value="")
        self.hv_notice_var = tk.StringVar(value="H/V is observational only and does not influence transition detection.")
        self._build()
        self.after(100, self._poll)

    def _build(self) -> None:
        workflow_tabs = ttk.Notebook(self)
        workflow_tabs.pack(fill="both", expand=True)
        single_pair_tab = ttk.Frame(workflow_tabs, padding=7)
        self.multi_pair_panel = MultiPairAnalysisPanel(
            workflow_tabs,
            get_genes=self.get_genes,
            get_quantiles=self.get_quantiles,
            get_scoring_options=self.get_scoring_options,
            get_output=lambda: self.output_var.get(),
        )
        workflow_tabs.add(single_pair_tab, text="Single Pair / Existing")
        workflow_tabs.add(self.multi_pair_panel, text="Multi-Pair Pre/Post")

        top = ttk.LabelFrame(single_pair_tab, text="Comparative Spatial Transition Analysis — v0.6-beta", padding=10)
        top.pack(fill="x")
        ttk.Label(
            top,
            text=(
                "Compares sample-level SpatialTX summaries. Delta = Target - Reference. "
                "No spot-wise subtraction or direct spatial registration is performed."
            ),
            foreground="#7c2d12",
            wraplength=920,
            justify="left",
        ).pack(anchor="w", fill="x")

        settings = ttk.Panedwindow(single_pair_tab, orient="horizontal")
        settings.pack(fill="x", pady=(8, 0))
        inputs = ttk.LabelFrame(settings, text="1–2  Comparison mode and input", padding=9)
        analysis = ttk.LabelFrame(settings, text="3  Reused analysis settings", padding=9)
        settings.add(inputs, weight=3)
        settings.add(analysis, weight=2)

        mode_row = ttk.Frame(inputs)
        mode_row.pack(fill="x")
        ttk.Label(mode_row, text="Mode").pack(side="left")
        self.mode_combo = ttk.Combobox(
            mode_row,
            textvariable=self.mode_var,
            values=("pairwise", "paired", "unpaired", "manifest_batch"),
            state="readonly",
            width=18,
        )
        self.mode_combo.pack(side="left", padx=(6, 8))
        self.mode_combo.bind("<<ComboboxSelected>>", lambda _event: self._update_mode_state())
        ttk.Button(mode_row, text="Refresh selected Main Mapper samples", command=self._refresh_samples).pack(side="left")

        pair = ttk.Frame(inputs)
        pair.pack(fill="x", pady=(7, 0))
        ttk.Label(pair, text="Sample A (Reference)").grid(row=0, column=0, sticky="w")
        ttk.Label(pair, text="Sample B (Target)").grid(row=0, column=1, sticky="w", padx=(7, 0))
        self.sample_a_combo = ttk.Combobox(pair, textvariable=self.sample_a_display_var, state="readonly")
        self.sample_b_combo = ttk.Combobox(pair, textvariable=self.sample_b_display_var, state="readonly")
        self.sample_a_combo.grid(row=1, column=0, sticky="ew")
        self.sample_b_combo.grid(row=1, column=1, sticky="ew", padx=(7, 0))
        self.sample_a_combo.bind("<<ComboboxSelected>>", lambda _event: self._sync_sample_path("a"))
        self.sample_b_combo.bind("<<ComboboxSelected>>", lambda _event: self._sync_sample_path("b"))
        self.sample_a_path_entry = ttk.Entry(pair, textvariable=self.sample_a_var, state="readonly")
        self.sample_b_path_entry = ttk.Entry(pair, textvariable=self.sample_b_var, state="readonly")
        self.sample_a_path_entry.grid(row=2, column=0, sticky="ew", pady=(3, 0))
        self.sample_b_path_entry.grid(row=2, column=1, sticky="ew", padx=(7, 0), pady=(3, 0))
        pair.columnconfigure(0, weight=1)
        pair.columnconfigure(1, weight=1)

        manifest = ttk.Frame(inputs)
        manifest.pack(fill="x", pady=(7, 0))
        ttk.Label(manifest, text="Manifest").pack(side="left")
        self.manifest_entry = ttk.Entry(manifest, textvariable=self.manifest_var)
        self.manifest_entry.pack(side="left", fill="x", expand=True, padx=6)
        self.manifest_button = ttk.Button(manifest, text="Load manifest…", command=self._browse_manifest)
        self.manifest_button.pack(side="left")

        group = ttk.Frame(inputs)
        group.pack(fill="x", pady=(7, 0))
        ttk.Label(group, text="Reference group").grid(row=0, column=0, sticky="w")
        ttk.Label(group, text="Target group").grid(row=0, column=1, sticky="w", padx=(7, 0))
        ttk.Label(group, text="Group column").grid(row=0, column=2, sticky="w", padx=(7, 0))
        ttk.Label(group, text="Pair ID column").grid(row=0, column=3, sticky="w", padx=(7, 0))
        ttk.Entry(group, textvariable=self.reference_var, width=14).grid(row=1, column=0, sticky="ew")
        ttk.Entry(group, textvariable=self.target_var, width=14).grid(row=1, column=1, sticky="ew", padx=(7, 0))
        ttk.Entry(group, width=12, state="readonly").grid(row=1, column=2, sticky="ew", padx=(7, 0))
        ttk.Entry(group, width=12, state="readonly").grid(row=1, column=3, sticky="ew", padx=(7, 0))
        group.winfo_children()[-2].configure(state="normal")
        group.winfo_children()[-2].insert(0, "group")
        group.winfo_children()[-2].configure(state="readonly")
        group.winfo_children()[-1].configure(state="normal")
        group.winfo_children()[-1].insert(0, "pair_id")
        group.winfo_children()[-1].configure(state="readonly")
        for column in range(4):
            group.columnconfigure(column, weight=1)

        ttk.Label(
            analysis,
            text="C/S programs, thresholds, and robustness settings are read from the current Main Mapper controls.",
            foreground="#4b5563",
            wraplength=400,
        ).pack(anchor="w", fill="x")
        graph = ttk.Frame(analysis)
        graph.pack(fill="x", pady=(7, 0))
        ttk.Label(graph, text="Context graph").grid(row=0, column=0, sticky="w")
        ttk.Label(graph, text="K").grid(row=0, column=1, sticky="w", padx=(6, 0))
        ttk.Label(graph, text="Radius (optional)").grid(row=0, column=2, sticky="w", padx=(6, 0))
        ttk.Combobox(graph, textvariable=self.graph_method_var, values=("knn", "radius", "lattice"), state="readonly", width=10).grid(row=1, column=0, sticky="ew")
        ttk.Spinbox(graph, from_=1, to=50, textvariable=self.graph_k_var, width=6).grid(row=1, column=1, sticky="ew", padx=(6, 0))
        ttk.Entry(graph, textvariable=self.graph_radius_var, width=10).grid(row=1, column=2, sticky="ew", padx=(6, 0))
        graph.columnconfigure(0, weight=1)
        graph.columnconfigure(2, weight=1)
        context = ttk.Frame(analysis)
        context.pack(fill="x", pady=(7, 0))
        ttk.Checkbutton(context, text="H_expr observational context", variable=self.h_var).grid(row=0, column=0, sticky="w")
        ttk.Checkbutton(context, text="V_expr observational context", variable=self.v_var).grid(row=1, column=0, sticky="w")
        ttk.Label(context, text="Context smoothing").grid(row=0, column=1, sticky="w", padx=(8, 0))
        ttk.Combobox(context, textvariable=self.context_smoothing_var, values=("none", "graph_mean"), state="readonly", width=12).grid(row=1, column=1, sticky="w", padx=(8, 0))
        stats = ttk.Frame(analysis)
        stats.pack(fill="x", pady=(7, 0))
        ttk.Label(stats, text="Statistical test").grid(row=0, column=0, sticky="w")
        ttk.Label(stats, text="FDR").grid(row=0, column=1, sticky="w", padx=(6, 0))
        ttk.Label(stats, text="Seed").grid(row=0, column=2, sticky="w", padx=(6, 0))
        ttk.Combobox(stats, textvariable=self.test_var, values=("auto", "wilcoxon", "mannwhitney", "paired_t", "welch_t"), state="readonly", width=16).grid(row=1, column=0, sticky="ew")
        ttk.Label(stats, text="Benjamini-Hochberg", foreground="#4b5563").grid(row=1, column=1, sticky="w", padx=(6, 0))
        ttk.Spinbox(stats, from_=0, to=2147483647, textvariable=self.seed_var, width=12).grid(row=1, column=2, sticky="ew", padx=(6, 0))

        actions = ttk.LabelFrame(single_pair_tab, text="4  Validate and run", padding=8)
        actions.pack(fill="x", pady=(8, 0))
        ttk.Label(actions, text="Output root").pack(side="left")
        ttk.Entry(actions, textvariable=self.output_var).pack(side="left", fill="x", expand=True, padx=6)
        ttk.Button(actions, text="Browse…", command=self._browse_output).pack(side="left")
        self.validate_button = ttk.Button(actions, text="Validate inputs", command=self._validate_async)
        self.validate_button.pack(side="left", padx=(7, 0))
        self.run_button = ttk.Button(actions, text="Run comparison", style="Primary.TButton", command=self._run_async)
        self.run_button.pack(side="left", padx=(5, 0))
        self.cancel_button = ttk.Button(actions, text="Cancel safely", command=self._cancel, state="disabled")
        self.cancel_button.pack(side="left", padx=(5, 0))

        result_tabs = ttk.Notebook(single_pair_tab)
        result_tabs.pack(fill="both", expand=True, pady=(8, 0))
        summary_tab = ttk.Frame(result_tabs, padding=6)
        delta_tab = ttk.Frame(result_tabs, padding=6)
        transition_tab = ttk.Frame(result_tabs, padding=6)
        figure_tab = ttk.Frame(result_tabs, padding=6)
        result_tabs.add(summary_tab, text="Summary & samples")
        result_tabs.add(delta_tab, text="Metric changes")
        result_tabs.add(transition_tab, text="Regimes & warnings")
        result_tabs.add(figure_tab, text="Figures")

        self.summary_text = tk.Text(summary_tab, height=6, wrap="word", state="disabled", background="#f8fafc")
        self.summary_text.pack(fill="x")
        self.scale_banner = ttk.Label(
            summary_tab,
            textvariable=self.scale_banner_var,
            foreground="#9a3412",
            background="#fff7ed",
            wraplength=920,
            justify="left",
            padding=6,
        )
        self.scale_banner.pack(fill="x", pady=(5, 0))
        ttk.Label(
            summary_tab,
            textvariable=self.hv_notice_var,
            foreground="#4b5563",
            wraplength=920,
            justify="left",
        ).pack(fill="x", pady=(4, 0))
        self.sample_tree = self._table(summary_tab)
        self.delta_tree = self._table(delta_tab)
        self.transition_tree = self._table(transition_tab, height=7)
        self.warning_text = tk.Text(transition_tab, height=7, wrap="word", state="disabled", background="#fff7ed")
        self.warning_text.pack(fill="both", expand=True, pady=(6, 0))

        figure_controls = ttk.Frame(figure_tab)
        figure_controls.pack(fill="x")
        self.figure_combo = ttk.Combobox(figure_controls, textvariable=self.figure_var, state="disabled")
        self.figure_combo.pack(side="left", fill="x", expand=True)
        self.figure_combo.bind("<<ComboboxSelected>>", lambda _event: self._show_figure())
        self.open_figure_button = ttk.Button(figure_controls, text="Open figure", command=self._open_figure, state="disabled")
        self.open_figure_button.pack(side="left", padx=(6, 0))
        self.open_results_button = ttk.Button(figure_controls, text="Open results", command=self._open_results, state="disabled")
        self.open_results_button.pack(side="left", padx=(5, 0))
        self.figure_canvas = tk.Canvas(figure_tab, background="#111827", highlightthickness=0, height=310)
        self.figure_canvas.pack(fill="both", expand=True, pady=(6, 0))
        self.figure_canvas.bind("<Configure>", self._schedule_figure_resize)

        footer = ttk.Frame(single_pair_tab)
        footer.pack(fill="x", pady=(7, 0))
        self.progress = ttk.Progressbar(footer, mode="indeterminate", length=160)
        self.progress.pack(side="left")
        ttk.Label(footer, textvariable=self.status_var, wraplength=780).pack(side="left", padx=8)
        self._refresh_samples()
        self._update_mode_state()

    @staticmethod
    def _table(parent, height: int = 11) -> ttk.Treeview:
        frame = ttk.Frame(parent)
        frame.pack(fill="both", expand=True, pady=(6, 0))
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
    def _set_text(widget: tk.Text, value: str) -> None:
        widget.configure(state="normal")
        widget.delete("1.0", "end")
        widget.insert("1.0", value)
        widget.configure(state="disabled")

    @staticmethod
    def _fill_tree(tree: ttk.Treeview, table: pd.DataFrame, preferred: list[str]) -> None:
        tree.delete(*tree.get_children())
        columns = [column for column in preferred if column in table]
        if not columns:
            columns = list(table.columns[:12])
        tree.configure(columns=columns)
        for column in columns:
            tree.heading(column, text=column)
            tree.column(column, width=min(220, max(90, len(column) * 9)), stretch=False)
        for index, row in table.head(2000).iterrows():
            values = []
            for column in columns:
                value = row[column]
                if pd.isna(value):
                    value = ""
                elif isinstance(value, float):
                    value = f"{value:.5g}"
                values.append(value)
            tree.insert("", "end", iid=f"row_{index}", values=values)

    def _refresh_samples(self) -> None:
        sample_paths = list(dict.fromkeys(Path(path) for path in self.get_samples()))
        labels, self._sample_path_by_label = _sample_display_labels(sample_paths)
        self._sample_label_by_path = {path: label for label, path in self._sample_path_by_label.items()}
        self.sample_a_combo.configure(values=labels)
        self.sample_b_combo.configure(values=labels)

        current_a = self.sample_a_var.get()
        current_b = self.sample_b_var.get()
        self.sample_a_display_var.set(self._sample_label_by_path.get(current_a, labels[0] if labels else ""))
        self.sample_b_display_var.set(
            self._sample_label_by_path.get(current_b, labels[1] if len(labels) > 1 else "")
        )
        self._sync_sample_path("a")
        self._sync_sample_path("b")

    def _sync_sample_path(self, slot: str) -> None:
        display_var = self.sample_a_display_var if slot == "a" else self.sample_b_display_var
        path_var = self.sample_a_var if slot == "a" else self.sample_b_var
        path_entry = self.sample_a_path_entry if slot == "a" else self.sample_b_path_entry
        path_var.set(self._sample_path_by_label.get(display_var.get(), ""))
        self.after_idle(lambda entry=path_entry: entry.xview_moveto(1.0))

    def _update_mode_state(self) -> None:
        pairwise = self.mode_var.get() == "pairwise"
        self.sample_a_combo.configure(state="readonly" if pairwise else "disabled")
        self.sample_b_combo.configure(state="readonly" if pairwise else "disabled")
        self.manifest_entry.configure(state="disabled" if pairwise else "normal")
        self.manifest_button.configure(state="disabled" if pairwise else "normal")

    def _browse_manifest(self) -> None:
        value = filedialog.askopenfilename(title="Comparative manifest", filetypes=[("CSV", "*.csv")], parent=self)
        if value:
            self.manifest_var.set(value)

    def _browse_output(self) -> None:
        value = filedialog.askdirectory(title="Comparative output root", parent=self)
        if value:
            self.output_var.set(value)

    def _records(self) -> list[SampleRecord]:
        mode = self.mode_var.get()
        if mode != "pairwise":
            if not self.manifest_var.get().strip():
                raise ValueError("Choose a comparative manifest CSV.")
            return load_comparative_manifest(self.manifest_var.get())
        self._sync_sample_path("a")
        self._sync_sample_path("b")
        a = Path(self.sample_a_var.get())
        b = Path(self.sample_b_var.get())
        if not a.is_file() or not b.is_file():
            raise ValueError("Select two existing H5AD samples from the Main Mapper sample list.")
        if a.resolve() == b.resolve():
            raise ValueError("Sample A and Sample B must be different files.")
        a_id, b_id = a.stem, b.stem
        if a_id == b_id:
            a_id, b_id = f"{a_id}_A", f"{b_id}_B"
        return [
            SampleRecord(a_id, a.resolve(), self.reference_var.get().strip()),
            SampleRecord(b_id, b.resolve(), self.target_var.get().strip()),
        ]

    def _config(self) -> ComparativeConfig:
        c_genes, s_genes = self.get_genes()
        c_q, s_q, g_q = self.get_quantiles()
        scoring = self.get_scoring_options()
        scoring_options = asdict(scoring) if hasattr(scoring, "__dataclass_fields__") else dict(scoring)
        graph = {"method": self.graph_method_var.get(), "k": int(self.graph_k_var.get()), "weighting": "binary"}
        radius = self.graph_radius_var.get().strip()
        if radius:
            graph["radius"] = float(radius)
        return ComparativeConfig(
            mode=self.mode_var.get(),
            reference=self.reference_var.get().strip(),
            target=self.target_var.get().strip(),
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
            statistical_test=self.test_var.get(),
            seed=int(self.seed_var.get()),
        )

    def _start_busy(self, message: str) -> None:
        self.busy = True
        self.cancel_event = Event()
        self.progress.start(12)
        self.status_var.set(message)
        self.validate_button.configure(state="disabled")
        self.run_button.configure(state="disabled")
        self.cancel_button.configure(state="normal")

    def _finish_busy(self) -> None:
        self.busy = False
        self.cancel_event = None
        self.progress.stop()
        self.validate_button.configure(state="normal")
        self.run_button.configure(state="normal")
        self.cancel_button.configure(state="disabled")

    def _validate_async(self) -> None:
        if self.busy:
            return
        try:
            records, config = self._records(), self._config()
            validate_record_structure(records, config)
        except Exception as exc:
            messagebox.showwarning("Comparative input validation", str(exc), parent=self)
            return
        self._start_busy("Validating H5AD files and spatial coordinates...")

        def work() -> None:
            try:
                table = preflight_records(records, config)
                self.events.put(("validated", table))
            except Exception as exc:
                self.events.put(("error", "Comparative validation failed", exc))

        threading.Thread(target=work, daemon=True).start()

    def _run_async(self) -> None:
        if self.busy:
            return
        try:
            records, config = self._records(), self._config()
            output = self.output_var.get().strip()
            if not output:
                raise ValueError("Choose a comparative output root.")
            validate_record_structure(records, config)
        except Exception as exc:
            messagebox.showwarning("Comparative Analysis", str(exc), parent=self)
            return
        self._start_busy("Running Comparative Spatial Transition Analysis...")

        def progress(message: str) -> None:
            self.events.put(("progress", message))

        def work() -> None:
            try:
                result = run_comparative_analysis(records, config, output, progress, self.cancel_event)
                self.events.put(("done", result))
            except Exception as exc:
                self.events.put(("error", "Comparative Analysis failed", exc))

        threading.Thread(target=work, daemon=True).start()

    def _cancel(self) -> None:
        if self.cancel_event is not None:
            self.cancel_event.set()
            self.status_var.set("Cancellation requested; current sample will finish before safe partial export.")

    def _poll(self) -> None:
        try:
            while True:
                event = self.events.get_nowait()
                if event[0] == "progress":
                    self.status_var.set(str(event[1]))
                elif event[0] == "validated":
                    table = event[1]
                    self._finish_busy()
                    self._fill_tree(self.sample_tree, table, ["sample_id", "group", "pair_id", "validation_status", "validation_error", "n_spots", "n_genes", "C_genes_present", "S_genes_present"])
                    ok = int(table["validation_status"].eq("ok").sum())
                    self.status_var.set(f"Input validation complete: {ok}/{len(table)} valid. Invalid rows remain visible and are never silently skipped.")
                elif event[0] == "done":
                    self._finish_busy()
                    self._load_result(event[1])
                elif event[0] == "error":
                    _, title, exc = event
                    self._finish_busy()
                    self.status_var.set(str(exc))
                    messagebox.showerror(title, str(exc), parent=self)
        except queue.Empty:
            pass
        if self.winfo_exists():
            self.after(100, self._poll)

    def _load_result(self, result: ComparativeRunResult) -> None:
        self.last_result = result
        self._set_text(self.summary_text, result.summary_text)
        self._fill_tree(self.sample_tree, result.sample_metrics, ["sample_id", "group", "QC_flag", "regime_label", "regime_confidence", "n_valid_spots", "tissue_component_count", "localized_interface_fraction", "diffuse_fraction", "transition_burden_score", "H_expr_available", "V_expr_available"])
        self._fill_tree(
            self.delta_tree,
            result.metric_change_table if not result.metric_change_table.empty else result.delta_metrics,
            [
                "comparison_id", "display_name", "reference_value", "target_value", "raw_delta",
                "standardized_delta", "normalized_delta", "symmetric_percent_change", "scale_sensitive",
                "observational_only", "warning",
            ],
        )
        self._fill_tree(self.transition_tree, result.regime_transitions, ["comparison_id", "reference_regime", "target_regime", "regime_transition", "transition_confidence_flag"])
        warning_lines = [f"- {row.get('sample_id', '')}: {row.get('message', '')}" for _, row in result.warnings.iterrows()]
        self._set_text(self.warning_text, "\n".join(warning_lines) if warning_lines else "No warnings recorded.")
        scale_cautions = (
            result.scale_warnings.loc[result.scale_warnings.get("severity", pd.Series(dtype=str)).eq("caution"), "message"]
            if not result.scale_warnings.empty and "message" in result.scale_warnings
            else pd.Series(dtype=str)
        )
        self.scale_banner_var.set(
            " ".join(dict.fromkeys(str(message) for message in scale_cautions if str(message).strip()))
            if len(scale_cautions)
            else "Sample-scale context: no configured substantial difference was detected."
        )
        self.hv_notice_var.set(
            "H/V is observational only. Centered H/V means are excluded from primary interpretation; "
            "available non-centered summaries use a pooled reference/target threshold."
        )
        figure_groups = {
            "comparative_summary_card.png": "Overview — Summary card",
            "comparative_group_distributions.png": "Overview — Group distributions",
            "comparative_metric_changes.png": "Overview — Raw delta (two independent panels)",
            "comparative_metric_changes_standardized.png": "Overview — Standardized change (group analyses)",
            "comparative_program_score_changes.png": "Program scores — Changes",
            "comparative_transition_changes.png": "Transition metrics — Changes",
            "comparative_graph_changes.png": "Graph metrics — Changes",
            "comparative_topology_raw_counts.png": "Topology: raw counts — Changes",
            "comparative_topology_normalized.png": "Topology: normalized — Changes",
            "comparative_sample_scale.png": "Sample scale — Reference vs target",
            "comparative_relative_changes.png": "Relative changes — Symmetric percent",
            "comparative_standardized_heatmap.png": "Standardized heatmap — Category grouped",
            "comparative_regime_transitions.png": "Regime transition — Descriptive",
            "comparative_side_by_side_maps.png": "Side-by-side maps — R = C - S",
            "comparative_HV_summary.png": "Observation layer: H/V — Summary",
            "comparative_H_V_context.png": "Observation layer: H/V — Compatibility output",
        }
        self.figure_label_to_path = {
            figure_groups.get(path.name, f"Other — {path.name}"): path for path in result.figures
        }
        labels = list(self.figure_label_to_path)
        self.figure_combo.configure(values=labels, state="readonly" if labels else "disabled")
        if labels:
            self.figure_var.set(labels[0])
            self._show_figure()
        self.open_figure_button.configure(state="normal" if labels else "disabled")
        self.open_results_button.configure(state="normal")
        ok = int(result.run_manifest["status"].eq("ok").sum())
        self.status_var.set(f"Comparative Analysis complete: {ok}/{len(result.run_manifest)} samples. Results: {result.run_dir}")

    def _current_figure(self) -> Path | None:
        if self.last_result is None:
            return None
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
        path = self._current_figure()
        self.figure_canvas.delete("all")
        if path is None or not path.is_file():
            self.figure_canvas.create_text(20, 20, anchor="nw", fill="#cbd5e1", text="Comparative figure unavailable.")
            self.current_photo = None
            return
        try:
            with Image.open(path) as source:
                source.load()
                display = source.convert("RGBA")
            display.thumbnail((max(100, self.figure_canvas.winfo_width() - 20), max(100, self.figure_canvas.winfo_height() - 20)), Image.Resampling.LANCZOS)
            self.current_photo = ImageTk.PhotoImage(display, master=self.figure_canvas)
            self.figure_canvas.create_image(max(1, self.figure_canvas.winfo_width()) // 2, max(1, self.figure_canvas.winfo_height()) // 2, image=self.current_photo)
        except (OSError, ValueError) as exc:
            self.current_photo = None
            self.figure_canvas.create_text(20, 20, anchor="nw", fill="#cbd5e1", text=f"Could not display figure: {exc}")

    @staticmethod
    def _open_path(path: Path) -> None:
        if os.name == "nt":
            os.startfile(path)  # type: ignore[attr-defined]
        else:
            subprocess.Popen(["open" if os.uname().sysname == "Darwin" else "xdg-open", str(path)])

    def _open_figure(self) -> None:
        path = self._current_figure()
        if path and path.is_file():
            try:
                self._open_path(path)
            except OSError as exc:
                messagebox.showerror("Open comparative figure", str(exc), parent=self)

    def _open_results(self) -> None:
        if self.last_result and self.last_result.run_dir.is_dir():
            try:
                self._open_path(self.last_result.run_dir)
            except OSError as exc:
                messagebox.showerror("Open comparative results", str(exc), parent=self)
