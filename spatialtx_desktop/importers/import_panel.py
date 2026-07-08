from __future__ import annotations

import os
import queue
import threading
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, ttk
from typing import Callable

from .mex_to_h5ad import convert_mex_to_h5ad, detect_mex_sample
from .validate_h5ad import validate_h5ad
from .visium_to_h5ad import convert_visium_to_h5ad, detect_visium_sample


class ImportConvertPanel(ttk.Frame):
    """Raw-format conversion boundary for the h5ad-centered Main Mapper."""

    def __init__(self, master, on_use_in_mapper: Callable[[Path], None]) -> None:
        super().__init__(master)
        ttk.Label(self, text="Import / Convert", font=("Segoe UI Semibold", 14)).pack(anchor="w")
        ttk.Label(
            self,
            text=(
                "Convert supported raw 10x inputs to H5AD first. The Main Mapper remains H5AD-centered "
                "and does not analyze raw H5, MEX/MTX, CSV, JSON, PNG, RDS, or parquet files directly."
            ),
            foreground="#4b5563",
            wraplength=780,
            justify="left",
        ).pack(anchor="w", fill="x", pady=(2, 10))

        tabs = ttk.Notebook(self)
        tabs.pack(fill="both", expand=True)
        mex_tab = ttk.Frame(tabs, padding=10)
        visium_tab = ttk.Frame(tabs, padding=10)
        tabs.add(mex_tab, text="Raw 10x MEX/MTX → H5AD")
        tabs.add(visium_tab, text="Raw Visium H5 + spatial → H5AD")
        ConverterSection(
            mex_tab,
            title="Raw 10x MEX/MTX → H5AD",
            description="matrix.mtx + barcodes.tsv + features.tsv/genes.tsv; .gz variants are supported.",
            detect=detect_mex_sample,
            convert=convert_mex_to_h5ad,
            require_spatial=False,
            on_use_in_mapper=on_use_in_mapper,
        ).pack(fill="both", expand=True)
        ConverterSection(
            visium_tab,
            title="Raw Visium H5 + spatial → H5AD",
            description=(
                "filtered_feature_bc_matrix.h5 plus tissue positions and scalefactors; optional tissue images "
                "and supported .gz variants are preserved."
            ),
            detect=detect_visium_sample,
            convert=convert_visium_to_h5ad,
            require_spatial=True,
            on_use_in_mapper=on_use_in_mapper,
        ).pack(fill="both", expand=True)


class ConverterSection(ttk.Frame):
    def __init__(
        self,
        master,
        *,
        title: str,
        description: str,
        detect,
        convert,
        require_spatial: bool,
        on_use_in_mapper: Callable[[Path], None],
    ) -> None:
        super().__init__(master)
        self.title = title
        self.detect_backend = detect
        self.convert_backend = convert
        self.require_spatial = require_spatial
        self.on_use_in_mapper = on_use_in_mapper
        self.last_output: Path | None = None
        self.busy = False
        self.messages: queue.Queue[tuple[str, object]] = queue.Queue()

        ttk.Label(self, text=title, font=("Segoe UI Semibold", 12)).pack(anchor="w")
        ttk.Label(self, text=description, foreground="#4b5563", wraplength=740, justify="left").pack(
            anchor="w", fill="x", pady=(2, 8)
        )

        form = ttk.LabelFrame(self, text="Input and canonical output", padding=10)
        form.pack(fill="x")
        self.input_var = tk.StringVar()
        self.output_var = tk.StringVar()
        self.name_var = tk.StringVar()
        self._folder_row(form, "Select input folder", self.input_var, self._browse_input)
        self._folder_row(form, "Select output folder", self.output_var, self._browse_output)
        row = ttk.Frame(form)
        row.pack(fill="x")
        ttk.Label(row, text="Sample name", width=20).pack(side="left")
        ttk.Entry(row, textvariable=self.name_var).pack(side="left", fill="x", expand=True)
        ttk.Label(row, text=".h5ad", foreground="#4b5563").pack(side="left", padx=(4, 0))

        actions = ttk.Frame(self)
        actions.pack(fill="x", pady=(9, 0))
        self.detect_button = ttk.Button(actions, text="Detect input", command=self._detect)
        self.detect_button.grid(row=0, column=0, sticky="w")
        self.convert_button = ttk.Button(actions, text="Convert to H5AD", command=self._convert)
        self.convert_button.grid(row=0, column=1, sticky="w", padx=(6, 0))
        self.validate_button = ttk.Button(actions, text="Validate converted H5AD", command=self._validate)
        self.validate_button.grid(row=0, column=2, sticky="w", padx=(6, 0))
        self.open_button = ttk.Button(actions, text="Open output folder", command=self._open_output_folder)
        self.open_button.grid(row=1, column=0, sticky="w", pady=(6, 0))
        self.mapper_button = ttk.Button(
            actions, text="Use in Main Mapper", command=self._use_in_mapper, state="disabled"
        )
        self.mapper_button.grid(row=1, column=1, sticky="w", padx=(6, 0), pady=(6, 0))
        self.progress = ttk.Progressbar(actions, mode="indeterminate", length=120)
        self.progress.grid(row=1, column=2, sticky="ew", padx=(8, 0), pady=(6, 0))
        actions.columnconfigure(2, weight=1)

        status_box = ttk.LabelFrame(self, text="Status log", padding=8)
        status_box.pack(fill="both", expand=True, pady=(9, 0))
        self.status_var = tk.StringVar(value="Select an input folder, output folder, and sample name.")
        ttk.Label(status_box, textvariable=self.status_var, wraplength=740, justify="left").pack(
            anchor="w", fill="x"
        )
        self.log = tk.Text(
            status_box,
            height=12,
            state="disabled",
            wrap="word",
            background="#f8fafc",
            padx=9,
            pady=7,
        )
        self.log.pack(fill="both", expand=True, pady=(6, 0))
        self.after(100, self._poll_messages)

    def _folder_row(self, parent, label: str, variable: tk.StringVar, command) -> None:
        row = ttk.Frame(parent)
        row.pack(fill="x", pady=(0, 7))
        ttk.Label(row, text=label, width=20).pack(side="left")
        ttk.Entry(row, textvariable=variable).pack(side="left", fill="x", expand=True)
        ttk.Button(row, text="Browse...", command=command).pack(side="left", padx=(6, 0))

    def _browse_input(self) -> None:
        value = filedialog.askdirectory(title=f"Select input folder for {self.title}", parent=self)
        if not value:
            return
        folder = Path(value)
        self.input_var.set(str(folder))
        self.name_var.set(folder.name)
        if not self.output_var.get().strip():
            self.output_var.set(str(folder.parent / "spatialtx_converted"))
        self._detect()

    def _browse_output(self) -> None:
        value = filedialog.askdirectory(title="Select H5AD output folder", parent=self)
        if value:
            self.output_var.set(value)

    def _write(self, text: str, *, clear: bool = False) -> None:
        self.log.configure(state="normal")
        if clear:
            self.log.delete("1.0", "end")
        if text:
            self.log.insert("end", text.rstrip() + "\n")
            self.log.see("end")
        self.log.configure(state="disabled")

    def _detect(self) -> dict | None:
        folder = self.input_var.get().strip()
        if not folder:
            messagebox.showwarning(self.title, "Select an input folder first.", parent=self)
            return None
        report = self.detect_backend(folder)
        self._write("", clear=True)
        if report["valid"]:
            suggested = report.get("suggested_sample_name")
            if suggested and self.name_var.get().strip() in {"", Path(folder).name}:
                self.name_var.set(str(suggested))
            self.status_var.set("Supported raw input detected. Ready to convert to H5AD.")
            for key, value in report.items():
                if key not in {"valid", "errors", "warnings", "inspection"} and value is not None:
                    self._write(f"{key}: {value}")
            if report.get("inspection"):
                for key, value in report["inspection"].items():
                    self._write(f"{key}: {value}")
            for warning in report.get("warnings", []):
                self._write(f"Warning: {warning}")
        else:
            self.status_var.set("The selected folder is not a supported complete input sample.")
            for error in report["errors"]:
                self._write(f"Error: {error}")
            for warning in report.get("warnings", []):
                self._write(f"Warning: {warning}")
        return report

    def _set_busy(self, value: bool) -> None:
        self.busy = value
        state = "disabled" if value else "normal"
        for button in (self.detect_button, self.convert_button, self.validate_button):
            button.configure(state=state)
        if value:
            self.mapper_button.configure(state="disabled")
            self.progress.start(12)
        else:
            self.progress.stop()

    def _convert(self) -> None:
        if self.busy:
            return
        report = self._detect()
        if report is None or not report["valid"]:
            messagebox.showerror(self.title, "Fix the input detection errors before conversion.", parent=self)
            return
        output_dir = self.output_var.get().strip()
        sample_name = self.name_var.get().strip()
        if not output_dir or not sample_name:
            messagebox.showwarning(self.title, "Select an output folder and sample name.", parent=self)
            return
        output = Path(output_dir).expanduser() / f"{sample_name}.h5ad"
        overwrite = False
        if output.exists():
            overwrite = messagebox.askyesno(
                self.title,
                f"This H5AD already exists:\n{output}\n\nReplace it?",
                parent=self,
            )
            if not overwrite:
                return
        self.last_output = None
        self._set_busy(True)
        self.status_var.set("Converting raw input to canonical H5AD...")

        def worker() -> None:
            try:
                result = self.convert_backend(
                    report.get("sample_folder", report.get("input_folder")),
                    output_dir,
                    sample_name,
                    overwrite=overwrite,
                    progress=lambda text: self.messages.put(("progress", text)),
                )
                self.messages.put(("success", result))
            except Exception as exc:
                self.messages.put(("error", exc))

        threading.Thread(target=worker, daemon=True).start()

    def _validate(self) -> None:
        initial = str(self.last_output) if self.last_output is not None else ""
        path = filedialog.askopenfilename(
            title="Validate converted H5AD",
            initialfile=Path(initial).name if initial else "",
            initialdir=str(Path(initial).parent) if initial else self.output_var.get(),
            filetypes=[("H5AD", "*.h5ad")],
            parent=self,
        )
        if not path:
            return
        report = validate_h5ad(path, require_spatial=self.require_spatial)
        self._write("Validation report:")
        for key, value in report.items():
            self._write(f"{key}: {value}")
        self.status_var.set("H5AD validation passed." if report["valid"] else "H5AD validation failed.")
        if report["valid"]:
            self.last_output = Path(path)
            self.mapper_button.configure(state="normal")
        else:
            messagebox.showerror(self.title, "H5AD validation failed. See the status log.", parent=self)

    def _open_output_folder(self) -> None:
        folder = self.last_output.parent if self.last_output is not None else Path(self.output_var.get()).expanduser()
        if not folder.is_dir():
            messagebox.showwarning(self.title, "The output folder does not exist yet.", parent=self)
            return
        os.startfile(folder)

    def _use_in_mapper(self) -> None:
        if self.last_output is None or not self.last_output.is_file():
            messagebox.showwarning(self.title, "Convert or validate an H5AD first.", parent=self)
            return
        self.on_use_in_mapper(self.last_output.parent)

    def _poll_messages(self) -> None:
        try:
            while True:
                kind, payload = self.messages.get_nowait()
                if kind == "progress":
                    self.status_var.set(str(payload))
                    self._write(str(payload))
                elif kind == "success":
                    path, report = payload
                    self.last_output = Path(path)
                    self._set_busy(False)
                    self.mapper_button.configure(state="normal")
                    self.status_var.set(f"Conversion complete and validated: {path}")
                    self._write(f"Success: {path}")
                    for key, value in report.items():
                        self._write(f"{key}: {value}")
                    messagebox.showinfo(
                        self.title,
                        f"Converted and validated:\n{path}\n\n"
                        f"{report['n_obs']} observations x {report['n_vars']} features",
                        parent=self,
                    )
                elif kind == "error":
                    self._set_busy(False)
                    self.status_var.set("Conversion failed. No completed H5AD was produced.")
                    self._write(f"Failure: {payload}")
                    messagebox.showerror(self.title, str(payload), parent=self)
        except queue.Empty:
            pass
        if self.winfo_exists():
            self.after(100, self._poll_messages)
