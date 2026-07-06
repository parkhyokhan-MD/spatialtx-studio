from __future__ import annotations

import os
import queue
import subprocess
import threading
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, ttk
from typing import Callable

import numpy as np
import pandas as pd

from .advanced_analysis import MODULE_LABELS, run_advanced_batch


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
        self.events: queue.Queue[tuple] = queue.Queue()
        self.busy = False
        self.last_run: Path | None = None
        self.dashboard_records: dict[str, dict] = {}
        self.selected_record: dict | None = None
        self.dashboard_counter = 0
        self._build()
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
        dashboard = ttk.Frame(self.tabs, padding=10)
        self.dashboard_tab = dashboard
        self.tabs.add(composition, text="Gene Composition")
        self.tabs.add(enrichment, text="Interface Enrichment")
        self.tabs.add(interaction, text="Cx/Sx Interaction")
        self.tabs.add(dashboard, text="Results Dashboard")
        self._build_composition(composition)
        self._build_enrichment(enrichment)
        self._build_interaction(interaction)
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
            "Calculates the relative contribution of every requested Cx and Sx gene on the same transformed expression scale used by v0.1. Missing genes remain visible and are marked explicitly.",
            "gene_composition.csv, 300-dpi PNG, vector PDF, analysis_metadata.json, and a run manifest.",
        )
        ttk.Button(parent, text="Run Gene Composition", style="Primary.TButton", command=lambda: self._run("composition")).pack(fill="x", pady=(14, 0))

    def _build_enrichment(self, parent) -> None:
        self._description(
            parent,
            "Interface-like versus non-interface composition",
            "Uses the unchanged v0.1 interface call. Reports group means, composition percentages, fold enrichment, Hedges' g, two-sided Mann-Whitney p-values, and Benjamini-Hochberg FDR.",
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
        ttk.Spinbox(settings, from_=0, to=9999, textvariable=self.permutations_var, width=10).grid(row=0, column=1, sticky="w", padx=(6, 18))
        ttk.Label(settings, text="Seed").grid(row=0, column=2, sticky="w")
        ttk.Spinbox(settings, from_=0, to=2147483647, textvariable=self.seed_var, width=14).grid(row=0, column=3, sticky="w", padx=6)
        ttk.Label(
            settings,
            text="Sx values are permuted across fixed coordinates, then local neighborhoods are recomputed.",
            foreground="#4b5563",
        ).grid(row=1, column=0, columnspan=4, sticky="w", pady=(7, 0))
        ttk.Button(parent, text="Run Cx/Sx Interaction", style="Primary.TButton", command=lambda: self._run("interaction")).pack(fill="x", pady=(14, 0))

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

    def _set_enabled(self, enabled: bool) -> None:
        for child in self.winfo_children():
            self._set_widget_enabled(child, enabled)

    def _set_widget_enabled(self, widget, enabled: bool) -> None:
        if isinstance(widget, (ttk.Button, ttk.Entry, ttk.Spinbox, ttk.Notebook)):
            widget.state(["!disabled"] if enabled else ["disabled"])
        for child in widget.winfo_children():
            self._set_widget_enabled(child, enabled)

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
