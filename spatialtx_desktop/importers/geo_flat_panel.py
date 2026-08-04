from __future__ import annotations

import os
import queue
import subprocess
import threading
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, ttk
from typing import Callable

from .geo_flat_visium import (
    GeoFlatSample,
    convert_geo_flat_directory,
    scan_geo_flat_directory,
    write_comparative_manifest,
)


class GeoFlatImportPanel(ttk.Frame):
    """User-selected, read-only-source workflow for flat standard Visium files."""

    def __init__(self, master, *, on_use_in_mapper: Callable[[Path], None], on_comparative_manifest=None) -> None:
        super().__init__(master)
        self.on_use_in_mapper = on_use_in_mapper
        self.on_comparative_manifest = on_comparative_manifest
        self.samples: list[GeoFlatSample] = []
        self.iid_to_sample: dict[str, GeoFlatSample] = {}
        self.converted_paths: dict[str, Path] = {}
        self.events: queue.Queue[tuple] = queue.Queue()
        self.busy = False
        self.source_var = tk.StringVar()
        self.output_var = tk.StringVar()
        self.recursive_var = tk.BooleanVar(value=False)
        self.status_var = tk.StringVar(value="Choose a source directory, then scan. Nothing is converted automatically.")
        self._build()
        self.after(100, self._poll)

    def _build(self) -> None:
        ttk.Label(self, text="GEO Flat Visium Directory", font=("Segoe UI Semibold", 12)).pack(anchor="w")
        ttk.Label(
            self,
            text=(
                "Groups supported files by exact full filename prefix. The selected source directory is read-only; "
                "H5AD files and reports are written only to the user-selected output directory."
            ),
            foreground="#4b5563", wraplength=880, justify="left",
        ).pack(anchor="w", fill="x", pady=(2, 8))

        source = ttk.LabelFrame(self, text="1  Source directory", padding=8)
        source.pack(fill="x")
        row = ttk.Frame(source); row.pack(fill="x")
        ttk.Entry(row, textvariable=self.source_var).pack(side="left", fill="x", expand=True)
        ttk.Button(row, text="Browse...", command=self._browse_source).pack(side="left", padx=(6, 0))
        self.scan_button = ttk.Button(row, text="Scan directory", command=self._scan)
        self.scan_button.pack(side="left", padx=(5, 0))
        ttk.Checkbutton(source, text="Include subdirectories (off by default)", variable=self.recursive_var).pack(
            anchor="w", pady=(5, 0)
        )

        box = ttk.LabelFrame(self, text="2  Detected inventory — select valid samples", padding=7)
        box.pack(fill="both", expand=True, pady=(8, 0))
        columns = (
            "selected", "prefix", "subject", "condition", "matrix", "barcodes", "features",
            "positions", "scalefactors", "images", "status", "warnings",
        )
        self.tree = ttk.Treeview(box, columns=columns, show="headings", selectmode="extended", height=10)
        labels = (
            "Select", "Sample prefix", "Parsed subject", "Parsed condition", "Matrix", "Barcodes",
            "Features", "Positions", "Scalefactors", "Images", "Status", "Warnings / errors",
        )
        widths = (55, 210, 105, 90, 115, 100, 100, 115, 115, 55, 75, 300)
        for column, label, width in zip(columns, labels, widths):
            self.tree.heading(column, text=label)
            self.tree.column(column, width=width, stretch=False)
        yscroll = ttk.Scrollbar(box, orient="vertical", command=self.tree.yview)
        xscroll = ttk.Scrollbar(box, orient="horizontal", command=self.tree.xview)
        self.tree.configure(yscrollcommand=yscroll.set, xscrollcommand=xscroll.set)
        self.tree.grid(row=0, column=0, sticky="nsew"); yscroll.grid(row=0, column=1, sticky="ns")
        xscroll.grid(row=1, column=0, sticky="ew")
        box.rowconfigure(0, weight=1); box.columnconfigure(0, weight=1)
        self.tree.bind("<<TreeviewSelect>>", lambda _event: self._mark_selection())
        select_row = ttk.Frame(box); select_row.grid(row=2, column=0, sticky="w", pady=(5, 0))
        ttk.Button(select_row, text="Select valid", command=self._select_valid).pack(side="left")
        ttk.Button(select_row, text="Clear", command=lambda: self.tree.selection_remove(self.tree.selection())).pack(
            side="left", padx=5
        )

        output = ttk.LabelFrame(self, text="3  Output and conversion", padding=8)
        output.pack(fill="x", pady=(8, 0))
        outrow = ttk.Frame(output); outrow.pack(fill="x")
        ttk.Label(outrow, text="Output folder", width=14).pack(side="left")
        ttk.Entry(outrow, textvariable=self.output_var).pack(side="left", fill="x", expand=True)
        ttk.Button(outrow, text="Browse...", command=self._browse_output).pack(side="left", padx=(6, 0))
        actions = ttk.Frame(output); actions.pack(fill="x", pady=(6, 0))
        self.convert_button = ttk.Button(actions, text="Convert selected to H5AD", command=self._convert)
        self.convert_button.pack(side="left")
        self.manifest_button = ttk.Button(actions, text="Review comparative manifest", command=self._review_manifest, state="disabled")
        self.manifest_button.pack(side="left", padx=5)
        self.mapper_button = ttk.Button(actions, text="Use in Main Mapper", command=self._use_in_mapper, state="disabled")
        self.mapper_button.pack(side="left")
        ttk.Button(actions, text="Open output", command=self._open_output).pack(side="left", padx=5)
        self.progress = ttk.Progressbar(actions, mode="indeterminate", length=130); self.progress.pack(side="right")

        status = ttk.LabelFrame(self, text="Progress and validation log", padding=7)
        status.pack(fill="both", expand=True, pady=(8, 0))
        ttk.Label(status, textvariable=self.status_var, wraplength=880, justify="left").pack(anchor="w", fill="x")
        self.log = tk.Text(status, height=7, state="disabled", wrap="word", background="#f8fafc")
        self.log.pack(fill="both", expand=True, pady=(5, 0))

    def _browse_source(self) -> None:
        value = filedialog.askdirectory(title="Select a GEO-style flat Visium directory", parent=self)
        if value:
            source = Path(value)
            self.source_var.set(str(source))
            if not self.output_var.get().strip():
                self.output_var.set(str(source.parent / f"{source.name}_spatialtx_converted"))

    def _browse_output(self) -> None:
        value = filedialog.askdirectory(title="Select GEO flat H5AD output directory", parent=self)
        if value:
            self.output_var.set(value)

    def _write(self, message: str, *, clear: bool = False) -> None:
        self.log.configure(state="normal")
        if clear:
            self.log.delete("1.0", "end")
        if message:
            self.log.insert("end", message.rstrip() + "\n"); self.log.see("end")
        self.log.configure(state="disabled")

    def _set_busy(self, busy: bool, message: str = "") -> None:
        self.busy = busy
        state = "disabled" if busy else "normal"
        self.scan_button.configure(state=state); self.convert_button.configure(state=state)
        self.progress.start(12) if busy else self.progress.stop()
        if message:
            self.status_var.set(message)

    def _scan(self) -> None:
        if self.busy:
            return
        source = self.source_var.get().strip()
        if not source:
            messagebox.showwarning("GEO Flat Visium Directory", "Select a source folder first.", parent=self); return
        recursive = bool(self.recursive_var.get())
        self._set_busy(True, "Scanning exact suffixes and validating sample prefixes..."); self._write("", clear=True)

        def worker() -> None:
            try:
                self.events.put(("scanned", scan_geo_flat_directory(source, recursive=recursive)))
            except Exception as exc:
                self.events.put(("error", "Directory scan failed", exc))

        threading.Thread(target=worker, daemon=True).start()

    def _load_inventory(self, samples: list[GeoFlatSample]) -> None:
        self.samples = samples; self.iid_to_sample.clear(); self.converted_paths.clear()
        self.tree.delete(*self.tree.get_children())
        self.manifest_button.configure(state="disabled"); self.mapper_button.configure(state="disabled")
        for index, sample in enumerate(samples):
            iid = f"sample_{index}"; self.iid_to_sample[iid] = sample
            messages = [*sample.warnings, *sample.errors]
            self.tree.insert("", "end", iid=iid, values=(
                "", sample.sample_prefix, sample.parsed_subject_id, sample.parsed_condition,
                sample.matrix_file.name if sample.matrix_file else "missing",
                sample.barcodes_file.name if sample.barcodes_file else "missing",
                sample.features_file.name if sample.features_file else "missing",
                sample.positions_file.name if sample.positions_file else "missing",
                sample.scalefactors_file.name if sample.scalefactors_file else "missing",
                len(sample.image_files), sample.validation_status, "; ".join(messages),
            ))
            self._write(f"{sample.sample_prefix}: {sample.validation_status}")
            for message in messages:
                self._write(f"  - {message}")
        valid = sum(sample.valid for sample in samples)
        self.status_var.set(
            f"Detected {len(samples)} prefix(es): {valid} convertible. Select explicitly; no pairing or conversion is automatic."
        )

    def _mark_selection(self) -> None:
        selected = set(self.tree.selection())
        for iid in self.tree.get_children():
            values = list(self.tree.item(iid, "values"))
            if values:
                values[0] = "Yes" if iid in selected else ""; self.tree.item(iid, values=values)

    def _select_valid(self) -> None:
        self.tree.selection_set([iid for iid, sample in self.iid_to_sample.items() if sample.valid])

    def _convert(self) -> None:
        if self.busy:
            return
        selected = [self.iid_to_sample[iid] for iid in self.tree.selection() if iid in self.iid_to_sample]
        if not selected:
            messagebox.showwarning("GEO Flat Visium Directory", "Select at least one valid sample.", parent=self); return
        invalid = [sample.sample_prefix for sample in selected if not sample.valid]
        if invalid:
            messagebox.showerror("GEO Flat Visium Directory", "Invalid samples: " + ", ".join(invalid), parent=self); return
        output = self.output_var.get().strip()
        if not output:
            messagebox.showwarning("GEO Flat Visium Directory", "Select an output folder.", parent=self); return
        existing = [Path(output) / f"{sample.sample_prefix}.h5ad" for sample in selected if (Path(output) / f"{sample.sample_prefix}.h5ad").exists()]
        overwrite = False
        if existing:
            overwrite = messagebox.askyesno(
                "Existing H5AD output",
                "These files exist:\n" + "\n".join(map(str, existing)) + "\n\nReplace these H5AD files?",
                parent=self,
            )
            if not overwrite:
                return
        source = self.source_var.get()
        recursive = bool(self.recursive_var.get())
        self._set_busy(True, f"Converting {len(selected)} selected sample(s)...")

        def worker() -> None:
            try:
                result = convert_geo_flat_directory(
                    source, output,
                    selected_samples=[sample.sample_prefix for sample in selected],
                    recursive=recursive, overwrite=overwrite,
                    progress=lambda message: self.events.put(("progress", message)),
                )
                self.events.put(("converted", result))
            except Exception as exc:
                self.events.put(("error", "GEO flat conversion failed", exc))

        threading.Thread(target=worker, daemon=True).start()

    def _review_manifest(self) -> None:
        if self.converted_paths:
            ManifestReviewDialog(
                self, [sample for sample in self.samples if sample.sample_prefix in self.converted_paths],
                self.converted_paths, Path(self.output_var.get()), self._manifest_saved,
            )

    def _manifest_saved(self, path: Path, confirmed: bool) -> None:
        self._write(f"Reviewed comparative manifest: {path}")
        if confirmed and self.on_comparative_manifest and messagebox.askyesno(
            "Comparative Analysis handoff", "Open this manifest in Comparative Analysis? Nothing will run automatically.", parent=self
        ):
            self.on_comparative_manifest(path)

    def _use_in_mapper(self) -> None:
        if self.converted_paths:
            self.on_use_in_mapper(Path(self.output_var.get()).expanduser().resolve())

    @staticmethod
    def _open_path(path: Path) -> None:
        if os.name == "nt":
            os.startfile(path)  # type: ignore[attr-defined]
        else:
            import platform
            subprocess.Popen(["open" if platform.system() == "Darwin" else "xdg-open", str(path)])

    def _open_output(self) -> None:
        path = Path(self.output_var.get()).expanduser()
        if not path.is_dir():
            messagebox.showwarning("GEO Flat Visium Directory", "The output folder does not exist.", parent=self); return
        try:
            self._open_path(path.resolve())
        except OSError as exc:
            messagebox.showerror("Open output", str(exc), parent=self)

    def _poll(self) -> None:
        try:
            while True:
                event = self.events.get_nowait()
                if event[0] == "progress":
                    self.status_var.set(str(event[1])); self._write(str(event[1]))
                elif event[0] == "scanned":
                    self._set_busy(False); self._load_inventory(event[1])
                elif event[0] == "converted":
                    result = event[1]; self._set_busy(False); self.converted_paths = result["converted_paths"]
                    summary = result["summary"]
                    for _, row in summary.iterrows():
                        self._write(f"{row['sample_prefix']}: {row['status']} {row.get('error_message', '')}")
                    success = len(self.converted_paths)
                    self.status_var.set(
                        f"Finished: {success}/{len(summary)} successful. Reports: {result['report_dir']}"
                    )
                    self.manifest_button.configure(state="normal" if success else "disabled")
                    self.mapper_button.configure(state="normal" if success else "disabled")
                    if len(summary) and not summary["status"].eq("success").all():
                        messagebox.showwarning("Conversion completed with failures", "Review the log and conversion summary.", parent=self)
                    elif success:
                        messagebox.showinfo("Conversion complete", f"Converted and validated {success} H5AD sample(s).", parent=self)
                elif event[0] == "error":
                    self._set_busy(False); _, title, exc = event
                    self.status_var.set(str(exc)); self._write(f"Error: {exc}"); messagebox.showerror(title, str(exc), parent=self)
        except queue.Empty:
            pass
        if self.winfo_exists():
            self.after(100, self._poll)


class ManifestReviewDialog(tk.Toplevel):
    """Editable approval gate for provisional filename-derived grouping."""

    def __init__(self, parent, samples, converted_paths, output_dir: Path, on_saved) -> None:
        super().__init__(parent)
        self.title("Review comparative manifest metadata"); self.geometry("1180x520"); self.transient(parent.winfo_toplevel())
        self.samples = samples; self.converted_paths = converted_paths; self.output_dir = output_dir; self.on_saved = on_saved
        self.entries: dict[str, dict[str, tk.StringVar]] = {}
        ttk.Label(
            self,
            text=(
                "Filename values are provisional. Review every group and pair ID before confirmation. "
                "Adjacent accessions are never paired automatically."
            ), foreground="#9a3412", wraplength=1120, justify="left",
        ).pack(anchor="w", fill="x", padx=12, pady=(12, 8))
        canvas = tk.Canvas(self, highlightthickness=0); scroll = ttk.Scrollbar(self, orient="vertical", command=canvas.yview)
        content = ttk.Frame(canvas); window = canvas.create_window((0, 0), window=content, anchor="nw")
        content.bind("<Configure>", lambda _event: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.bind("<Configure>", lambda event: canvas.itemconfigure(window, width=event.width))
        canvas.configure(yscrollcommand=scroll.set); canvas.pack(side="left", fill="both", expand=True, padx=(12, 0)); scroll.pack(side="right", fill="y", padx=(0, 12))
        headers = ("Sample ID", "Accession", "Parsed subject", "Group", "Pair ID", "Condition", "Batch", "Notes")
        for column, label in enumerate(headers):
            ttk.Label(content, text=label, font=("Segoe UI Semibold", 9)).grid(row=0, column=column, sticky="w", padx=3)
        for row, sample in enumerate(samples, 1):
            ttk.Label(content, text=sample.sample_prefix, wraplength=190).grid(row=row, column=0, sticky="w", padx=3, pady=3)
            ttk.Label(content, text=sample.geo_accession).grid(row=row, column=1, sticky="w", padx=3)
            ttk.Label(content, text=sample.parsed_subject_id).grid(row=row, column=2, sticky="w", padx=3)
            variables = {
                "group": tk.StringVar(value=sample.parsed_condition), "pair_id": tk.StringVar(value=sample.parsed_subject_id),
                "condition": tk.StringVar(value=sample.parsed_condition), "batch": tk.StringVar(),
                "notes": tk.StringVar(value="Reviewed filename-derived metadata."),
            }
            self.entries[sample.sample_prefix] = variables
            for column, field in enumerate(("group", "pair_id", "condition", "batch", "notes"), 3):
                ttk.Entry(content, textvariable=variables[field], width=30 if field == "notes" else 13).grid(row=row, column=column, sticky="ew", padx=3, pady=3)
        buttons = ttk.Frame(self); buttons.pack(fill="x", padx=12, pady=10)
        ttk.Button(buttons, text="Save unconfirmed draft", command=lambda: self._save(False)).pack(side="left")
        ttk.Button(buttons, text="Confirm reviewed metadata and save", command=lambda: self._save(True)).pack(side="left", padx=6)
        ttk.Button(buttons, text="Cancel", command=self.destroy).pack(side="right"); self.grab_set()

    def _save(self, confirmed: bool) -> None:
        if confirmed and not messagebox.askyesno(
            "Confirm filename-derived metadata", "Have you reviewed every group and pair ID?", parent=self
        ):
            return
        value = filedialog.asksaveasfilename(
            title="Save comparative manifest", initialdir=str(self.output_dir),
            initialfile="comparative_manifest_confirmed.csv" if confirmed else "comparative_manifest_draft.csv",
            defaultextension=".csv", filetypes=[("CSV", "*.csv")], parent=self,
        )
        if not value:
            return
        edits = {sample_id: {field: var.get() for field, var in values.items()} for sample_id, values in self.entries.items()}
        try:
            path = write_comparative_manifest(
                self.samples, self.converted_paths, value, confirmed=confirmed, edits=edits,
            )
        except Exception as exc:
            messagebox.showerror("Save comparative manifest", str(exc), parent=self); return
        self.grab_release(); self.destroy(); self.on_saved(path, confirmed)
