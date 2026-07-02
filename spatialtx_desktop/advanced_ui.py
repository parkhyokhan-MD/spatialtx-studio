from __future__ import annotations

import queue
import threading
import tkinter as tk
from datetime import datetime
from pathlib import Path
from tkinter import filedialog, messagebox, ttk

from .advanced import (
    annotate_sequences,
    build_ligrec_skeleton,
    build_read_evidence_plan,
    compare_pre_post,
    convert_mex,
    export_fasta_template,
    export_qubo_pool,
    filter_receptor_membrane,
    find_mex_folders,
    inspect_mex,
    scan_pre_post_pairs,
    validate_h5ad,
)


class AdvancedToolsPanel(ttk.Frame):
    def __init__(self, parent, on_qubo_pool=None) -> None:
        super().__init__(parent)
        self.on_qubo_pool = on_qubo_pool
        self.mex_folders: list[Path] = []
        self.latest_table: Path | None = None
        self.raw_pre = tk.StringVar()
        self.raw_post = tk.StringVar()
        self.pre_h5ad = tk.StringVar()
        self.post_h5ad = tk.StringVar()
        self.input_table = tk.StringVar()
        self.enabled_var = tk.BooleanVar(value=False)
        self.busy = False
        self.events: queue.Queue[tuple] = queue.Queue()
        self.output_var = tk.StringVar(value=str(Path.cwd() / "advanced_results"))
        self._build()
        self._toggle_enabled()
        self.after(100, self._poll_events)

    def _build(self) -> None:
        header = ttk.LabelFrame(self, text="0  Enable Advanced tools", padding=10)
        header.pack(fill="x")
        ttk.Checkbutton(header, text="I understand these tools are experimental", variable=self.enabled_var, command=self._toggle_enabled).pack(side="left")
        ttk.Label(header, text="Output folder").pack(side="left", padx=(18, 5))
        ttk.Entry(header, textvariable=self.output_var).pack(side="left", fill="x", expand=True)
        ttk.Button(header, text="Browse...", command=self._browse_output).pack(side="left", padx=(5, 0))

        self.body = ttk.Frame(self)
        self.body.pack(fill="both", expand=True, pady=(8, 0))
        self.tabs = ttk.Notebook(self.body)
        self.tabs.pack(fill="both", expand=True)
        raw = ttk.Frame(self.tabs, padding=10)
        candidates = ttk.Frame(self.tabs, padding=10)
        utilities = ttk.Frame(self.tabs, padding=10)
        self.tabs.add(raw, text="A  Raw MEX to h5ad")
        self.tabs.add(candidates, text="B  Pre/Post comparison")
        self.tabs.add(utilities, text="C  Ligand/Receptor utilities")
        self._build_raw(raw)
        self._build_candidates(candidates)
        self._build_utilities(utilities)

        footer = ttk.Frame(self)
        footer.pack(fill="x", pady=(8, 0))
        self.progress = ttk.Progressbar(footer, mode="indeterminate", length=160)
        self.progress.pack(side="left")
        self.status = ttk.Label(footer, text="Advanced tools are disabled.", foreground="#7c2d12")
        self.status.pack(side="left", padx=10)

    def _build_raw(self, parent) -> None:
        scan = ttk.LabelFrame(parent, text="A0B  Raw folder scan", padding=8)
        scan.pack(fill="x")
        self.raw_root = tk.StringVar()
        ttk.Entry(scan, textvariable=self.raw_root).pack(side="left", fill="x", expand=True)
        ttk.Button(scan, text="Root...", command=self._browse_raw_root).pack(side="left", padx=5)
        ttk.Button(scan, text="Scan MEX folders", command=self._scan_mex).pack(side="left")
        self.raw_tree = ttk.Treeview(parent, columns=("folder",), show="headings", height=7, selectmode="browse")
        self.raw_tree.heading("folder", text="Complete 10x MEX folders")
        self.raw_tree.column("folder", width=700)
        self.raw_tree.pack(fill="both", expand=True, pady=(8, 0))
        setrow = ttk.Frame(parent); setrow.pack(fill="x", pady=(6, 0))
        ttk.Button(setrow, text="Set selected folder -> Raw pre", command=lambda: self._set_raw("pre")).pack(side="left")
        ttk.Button(setrow, text="Set selected folder -> Raw post", command=lambda: self._set_raw("post")).pack(side="left", padx=5)
        ttk.Button(setrow, text="A0  Inspect selected", command=self._inspect_selected).pack(side="left")
        paths = ttk.LabelFrame(parent, text="Selected raw folders", padding=8); paths.pack(fill="x", pady=(8, 0))
        ttk.Label(paths, text="Pre").grid(row=0, column=0, sticky="w"); ttk.Entry(paths, textvariable=self.raw_pre).grid(row=0, column=1, sticky="ew", padx=5)
        ttk.Label(paths, text="Post").grid(row=1, column=0, sticky="w"); ttk.Entry(paths, textvariable=self.raw_post).grid(row=1, column=1, sticky="ew", padx=5, pady=(4, 0))
        paths.columnconfigure(1, weight=1)
        actions = ttk.Frame(parent); actions.pack(fill="x", pady=(8, 0))
        ttk.Button(actions, text="A2  Convert pre", command=lambda: self._convert("pre")).pack(side="left")
        ttk.Button(actions, text="A2  Convert post", command=lambda: self._convert("post")).pack(side="left", padx=5)
        ttk.Button(actions, text="A2B  Validate h5ad...", command=self._validate_file).pack(side="left")
        self.raw_report = tk.Text(parent, height=7, wrap="word", state="disabled", background="#f8fafc")
        self.raw_report.pack(fill="x", pady=(8, 0))

    def _build_candidates(self, parent) -> None:
        guide = ttk.LabelFrame(parent, text="A3 -> A4 -> A5  Optional hypothesis-generation workflow", padding=9)
        guide.pack(fill="x")
        ttk.Label(
            guide,
            text=(
                "A3 — Exploratory candidate comparison: compares two conditions (for example pre/post, control/treated, or region A/B) using "
                "normalized mean-expression contrast and detection-fraction change. Output genes are condition-associated exploratory candidates.\n\n"
                "A4 — Lightweight filtering: prioritizes receptor-like, membrane-associated, transporter-like, and surface-like candidates "
                "using gene-symbol heuristics. It is a literature-review and experiment-prioritization aid, not protein-function validation.\n\n"
                "A5 — QUBO handoff: saves the filtered table as a bounded candidate pool and loads its gene list into the C/S QUBO selector. "
                "QUBO then evaluates which k-gene combination best explains the requested spatial side."
            ),
            wraplength=780, justify="left",
        ).pack(anchor="w", fill="x")
        ttk.Label(
            guide,
            text=(
                "These tools do not validate drug response, receptor function, ligand-receptor binding, read-level evidence, "
                "or clinical biomarkers. A3 currently performs expression/detection contrast; spatial direction is evaluated downstream by QUBO."
            ),
            foreground="#9a3412", wraplength=780, justify="left",
        ).pack(anchor="w", fill="x", pady=(7, 0))
        pairscan = ttk.LabelFrame(parent, text="A1  Scan pre/post pairs", padding=8); pairscan.pack(fill="x")
        self.pair_root = tk.StringVar(value=str(Path.cwd()))
        ttk.Entry(pairscan, textvariable=self.pair_root).pack(side="left", fill="x", expand=True)
        ttk.Button(pairscan, text="Folder...", command=self._browse_pair_root).pack(side="left", padx=5)
        ttk.Button(pairscan, text="Scan pairs", command=self._scan_pairs).pack(side="left")
        self.pair_tree = ttk.Treeview(parent, columns=("pair", "pre", "post"), show="headings", height=5, selectmode="browse")
        for column, title, width in (("pair", "Pair", 120), ("pre", "Pre h5ad", 300), ("post", "Post h5ad", 300)):
            self.pair_tree.heading(column, text=title); self.pair_tree.column(column, width=width)
        self.pair_tree.pack(fill="x", pady=(8, 0)); self.pair_tree.bind("<<TreeviewSelect>>", self._use_pair)
        paths = ttk.LabelFrame(parent, text="Comparison inputs", padding=8); paths.pack(fill="x", pady=(8, 0))
        ttk.Label(paths, text="Pre h5ad").grid(row=0, column=0, sticky="w"); ttk.Entry(paths, textvariable=self.pre_h5ad).grid(row=0, column=1, sticky="ew", padx=5)
        ttk.Label(paths, text="Post h5ad").grid(row=1, column=0, sticky="w"); ttk.Entry(paths, textvariable=self.post_h5ad).grid(row=1, column=1, sticky="ew", padx=5, pady=(4, 0))
        paths.columnconfigure(1, weight=1)
        buttons = ttk.Frame(parent); buttons.pack(fill="x", pady=(8, 0))
        ttk.Button(buttons, text="A3  Candidate comparison", command=self._compare).pack(side="left")
        ttk.Button(buttons, text="A4  Receptor/membrane filter", command=self._filter).pack(side="left", padx=5)
        ttk.Button(buttons, text="A5  Export to QUBO", command=self._export_qubo).pack(side="left")
        self.candidate_report = tk.Text(parent, height=9, wrap="word", state="disabled", background="#f8fafc")
        self.candidate_report.pack(fill="both", expand=True, pady=(8, 0))

    def _build_utilities(self, parent) -> None:
        source = ttk.LabelFrame(parent, text="Candidate/annotation table", padding=8); source.pack(fill="x")
        ttk.Entry(source, textvariable=self.input_table).pack(side="left", fill="x", expand=True)
        ttk.Button(source, text="CSV...", command=self._browse_table).pack(side="left", padx=(5, 0))
        actions = ttk.LabelFrame(parent, text="Experimental utilities", padding=10); actions.pack(fill="x", pady=(8, 0))
        ttk.Button(actions, text="A6  Sequence annotation template", command=self._annotate).grid(row=0, column=0, sticky="ew")
        ttk.Button(actions, text="A7  Ligand/receptor skeleton", command=self._ligrec).grid(row=0, column=1, sticky="ew", padx=(6, 0))
        ttk.Button(actions, text="A8  FASTA/template export", command=self._fasta).grid(row=1, column=0, sticky="ew", pady=(6, 0))
        ttk.Button(actions, text="A9  Read-evidence plan", command=self._evidence).grid(row=1, column=1, sticky="ew", padx=(6, 0), pady=(6, 0))
        actions.columnconfigure(0, weight=1); actions.columnconfigure(1, weight=1)
        ttk.Label(
            parent,
            text="These utilities create local review templates and candidate skeletons. They do not query external sequence or interaction databases.",
            foreground="#7c2d12", wraplength=700,
        ).pack(anchor="w", pady=(8, 0))
        self.utility_report = tk.Text(parent, height=15, wrap="word", state="disabled", background="#f8fafc")
        self.utility_report.pack(fill="both", expand=True, pady=(8, 0))

    def _toggle_enabled(self) -> None:
        if self.busy:
            self.enabled_var.set(True)
            return
        enabled = self.enabled_var.get()
        self._set_descendant_state(self.body, enabled)
        self.status.configure(text="Advanced tools enabled." if enabled else "Advanced tools are disabled.", foreground="#166534" if enabled else "#7c2d12")

    def _set_descendant_state(self, widget, enabled: bool) -> None:
        for child in widget.winfo_children():
            if isinstance(child, (ttk.Button, ttk.Entry, ttk.Treeview, ttk.Spinbox, ttk.Combobox, ttk.Notebook)):
                child.state(["!disabled"] if enabled else ["disabled"])
            self._set_descendant_state(child, enabled)

    def _out(self) -> Path:
        path = Path(self.output_var.get()).expanduser().resolve(); path.mkdir(parents=True, exist_ok=True); return path

    def _new_output(self, filename: str) -> Path:
        path = self._out() / filename
        if not path.exists():
            return path
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
        return path.with_name(f"{path.stem}_{stamp}{path.suffix}")

    def _background(self, label: str, function, callback=None) -> None:
        if self.busy:
            messagebox.showinfo("Advanced tools", "Wait for the current Advanced job to finish.", parent=self)
            return
        self.busy = True
        self._set_descendant_state(self.body, False)
        self.progress.start(12); self.status.configure(text=label, foreground="#334155")
        def work():
            try:
                result = function()
                self.events.put(("done", label, result, callback))
            except Exception as exc:
                self.events.put(("error", label, exc))
        threading.Thread(target=work, daemon=True).start()

    def _poll_events(self) -> None:
        try:
            while True:
                event = self.events.get_nowait()
                if event[0] == "done":
                    _, label, result, callback = event
                    self._finished(label, result, callback)
                else:
                    _, label, exc = event
                    self._failed(label, exc)
        except queue.Empty:
            pass
        if self.winfo_exists():
            self.after(100, self._poll_events)

    def _finished(self, label, result, callback) -> None:
        self.busy = False
        self._set_descendant_state(self.body, self.enabled_var.get())
        self.progress.stop(); self.status.configure(text=f"Completed: {label}", foreground="#166534")
        if callback:
            try:
                callback(result)
            except Exception as exc:
                self._failed(f"{label} result handling", exc)

    def _failed(self, label, exc) -> None:
        self.busy = False
        self._set_descendant_state(self.body, self.enabled_var.get())
        self.progress.stop(); self.status.configure(text=f"Failed: {label}", foreground="#b91c1c")
        messagebox.showerror("Advanced tool failed", str(exc), parent=self)

    @staticmethod
    def _write(text_widget: tk.Text, value: str) -> None:
        text_widget.configure(state="normal"); text_widget.delete("1.0", "end"); text_widget.insert("1.0", value); text_widget.configure(state="disabled")

    def _browse_output(self):
        value = filedialog.askdirectory(title="Advanced output folder"); self.output_var.set(value or self.output_var.get())
    def _browse_raw_root(self):
        value = filedialog.askdirectory(title="Raw 10x root folder"); self.raw_root.set(value or self.raw_root.get())
    def _browse_pair_root(self):
        value = filedialog.askdirectory(title="Folder containing pre/post h5ad files"); self.pair_root.set(value or self.pair_root.get())
    def _browse_table(self):
        value = filedialog.askopenfilename(title="Select candidate CSV", filetypes=[("CSV", "*.csv")]); self.input_table.set(value or self.input_table.get())

    def _scan_mex(self):
        root = self.raw_root.get()
        self._background("Raw folder scan", lambda: find_mex_folders(root), self._show_mex)
    def _show_mex(self, folders):
        self.mex_folders = folders; self.raw_tree.delete(*self.raw_tree.get_children())
        for i, folder in enumerate(folders): self.raw_tree.insert("", "end", iid=str(i), values=(str(folder),))
        self._write(self.raw_report, f"Found {len(folders)} complete MEX folder(s).")
    def _selected_mex(self) -> Path:
        selected = self.raw_tree.selection()
        if not selected: raise ValueError("Select a scanned MEX folder first.")
        return self.mex_folders[int(selected[0])]
    def _set_raw(self, stage):
        try: (self.raw_pre if stage == "pre" else self.raw_post).set(str(self._selected_mex()))
        except Exception as exc: messagebox.showwarning("Raw folder", str(exc), parent=self)
    def _inspect_selected(self):
        try: folder = self._selected_mex()
        except Exception as exc: messagebox.showwarning("Inspect", str(exc), parent=self); return
        self._background("MEX inspection", lambda: inspect_mex(folder), lambda report: self._write(self.raw_report, "\n".join(f"{k}: {v}" for k, v in report.items())))
    def _convert(self, stage):
        source = self.raw_pre.get() if stage == "pre" else self.raw_post.get()
        if not source: messagebox.showwarning("Convert", f"Set the raw {stage} folder first.", parent=self); return
        output = self._new_output(f"{Path(source).name}_{stage}.h5ad")
        def done(result):
            path, report = result; (self.pre_h5ad if stage == "pre" else self.post_h5ad).set(str(path)); self.pair_root.set(str(path.parent))
            self._write(self.raw_report, "\n".join(f"{k}: {v}" for k, v in report.items()))
        self._background(f"Convert raw {stage}", lambda: convert_mex(source, output), done)
    def _validate_file(self):
        path = filedialog.askopenfilename(title="Validate h5ad", filetypes=[("h5ad", "*.h5ad")])
        if path: self._background("h5ad validation", lambda: validate_h5ad(path), lambda report: self._write(self.raw_report, "\n".join(f"{k}: {v}" for k, v in report.items())))

    def _scan_pairs(self):
        root = self.pair_root.get()
        self._background("Pre/post pair scan", lambda: scan_pre_post_pairs(root), self._show_pairs)
    def _show_pairs(self, table):
        self.pair_tree.delete(*self.pair_tree.get_children())
        for i, row in table.iterrows(): self.pair_tree.insert("", "end", iid=str(i), values=(row["pair"], row["pre_h5ad"], row["post_h5ad"]))
        self._write(self.candidate_report, f"Found {len(table)} complete pre/post pair(s).")
    def _use_pair(self, _event=None):
        selected = self.pair_tree.selection()
        if selected:
            values = self.pair_tree.item(selected[0], "values"); self.pre_h5ad.set(values[1]); self.post_h5ad.set(values[2])
    def _compare(self):
        pre_path, post_path = self.pre_h5ad.get(), self.post_h5ad.get()
        output = self._new_output("pre_post_candidate_comparison.csv")
        def done(table): self.latest_table = output; self.input_table.set(str(output)); self._write(self.candidate_report, f"Saved {len(table)} gene candidates:\n{output}\n\n{table.head(12).to_string(index=False)}")
        self._background("Pre/post candidate comparison", lambda: compare_pre_post(pre_path, post_path, output), done)
    def _filter(self):
        source = self.latest_table or Path(self.input_table.get()); output = self._new_output("receptor_membrane_candidates.csv")
        def done(table): self.latest_table = output; self.input_table.set(str(output)); self._write(self.candidate_report, f"Saved {len(table)} filtered candidates:\n{output}\n\n{table.head(12).to_string(index=False)}")
        self._background("Receptor/membrane filter", lambda: filter_receptor_membrane(source, output), done)
    def _export_qubo(self):
        source = self.latest_table or Path(self.input_table.get()); output = self._new_output("qubo_candidate_pool.csv")
        def done(genes):
            if self.on_qubo_pool: self.on_qubo_pool(genes, output)
            self._write(self.candidate_report, f"Exported {len(genes)} genes to the QUBO candidate pool:\n{output}\n\n" + ", ".join(genes))
        self._background("QUBO pool export", lambda: export_qubo_pool(source, output), done)

    def _source_table(self) -> Path:
        path = Path(self.input_table.get())
        if not path.is_file(): raise ValueError("Select a valid input CSV first.")
        return path
    def _annotate(self):
        try: source = self._source_table()
        except Exception as exc: messagebox.showwarning("Sequence annotation", str(exc), parent=self); return
        output = self._new_output("sequence_annotation_template.csv")
        self._background("Sequence annotation template", lambda: annotate_sequences(source, output), lambda table: self._utility_done(output, table))
    def _ligrec(self):
        try: source = self._source_table()
        except Exception as exc: messagebox.showwarning("Ligand/receptor skeleton", str(exc), parent=self); return
        output = self._new_output("ligand_receptor_skeleton.csv")
        self._background("Ligand/receptor skeleton", lambda: build_ligrec_skeleton(source, output), lambda table: self._utility_done(output, table))
    def _fasta(self):
        try: source = self._source_table()
        except Exception as exc: messagebox.showwarning("FASTA/template export", str(exc), parent=self); return
        prefix = self._out() / "spatialtx_sequences"
        if prefix.with_suffix(".fasta").exists() or prefix.with_suffix(".sequence_template.csv").exists():
            prefix = self._out() / f"spatialtx_sequences_{datetime.now().strftime('%Y%m%d_%H%M%S_%f')}"
        self._background("FASTA/template export", lambda: export_fasta_template(source, prefix), lambda result: self._write(self.utility_report, f"FASTA: {result[0]}\nTemplate: {result[1]}\nSequence records: {result[2]}"))
    def _evidence(self):
        try: source = self._source_table()
        except Exception as exc: messagebox.showwarning("Read-evidence plan", str(exc), parent=self); return
        output = self._new_output("read_evidence_plan.csv")
        self._background("Read-evidence plan", lambda: build_read_evidence_plan(source, output), lambda table: self._utility_done(output, table))
    def _utility_done(self, output, table):
        self.input_table.set(str(output)); self._write(self.utility_report, f"Saved {len(table)} row(s):\n{output}\n\n{table.head(12).to_string(index=False)}")
