from __future__ import annotations

import ast
import json
import math
import tkinter as tk
from dataclasses import dataclass
from pathlib import Path
from tkinter import messagebox, ttk
from typing import Callable, Iterable

import pandas as pd
from PIL import Image, ImageTk


INITIAL_MESSAGE = "Run Spatial Graph & Neighborhood to view generated context and graph figures."

# Long suffixes intentionally precede their shorter counterparts.
FIGURE_SUFFIX_LABELS = [
    ("_H_expr_unsmoothed_map.png", "H_expr unsmoothed"),
    ("_H_expr_smoothed_map.png", "H_expr graph-smoothed"),
    ("_V_expr_unsmoothed_map.png", "V_expr unsmoothed"),
    ("_V_expr_smoothed_map.png", "V_expr graph-smoothed"),
    ("_H_V_association_map.png", "H/V joint high-state map"),
    ("_neighborhood_enrichment_heatmap.png", "Neighborhood enrichment heatmap"),
    ("_H_expr_map.png", "H_expr hypoxia-associated expression field"),
    ("_V_expr_map.png", "V_expr endothelial/angiogenic expression proxy"),
    ("_graph_qc.png", "Graph QC"),
]

FIGURE_PRIORITY = {
    "_H_V_association_map.png": 0,
    "_H_expr_map.png": 1,
    "_V_expr_map.png": 2,
    "_neighborhood_enrichment_heatmap.png": 3,
    "_graph_qc.png": 4,
    "_H_expr_smoothed_map.png": 5,
    "_H_expr_unsmoothed_map.png": 6,
    "_V_expr_smoothed_map.png": 7,
    "_V_expr_unsmoothed_map.png": 8,
}


@dataclass(frozen=True)
class FigureRecord:
    sample: str
    path: Path
    label: str
    suffix: str
    priority: int

    @property
    def selector_label(self) -> str:
        return f"{self.label} — {self.path.name}"


@dataclass(frozen=True)
class ManifestSummary:
    successful_samples: tuple[str, ...]
    failed_samples: tuple[tuple[str, str], ...]
    total_rows: int
    warning: str = ""


@dataclass(frozen=True)
class ContextFieldRecord:
    field: str
    status: str
    matched_gene_count: str
    requested_gene_count: str
    coverage_fraction: str
    score_method: str
    smoothing_method: str
    high_state_fraction: str
    warnings: tuple[str, ...]


@dataclass(frozen=True)
class ContextSummary:
    available: bool
    records: tuple[ContextFieldRecord, ...]
    message: str = ""


def _is_missing(value) -> bool:
    if value is None:
        return True
    if isinstance(value, str):
        return not value.strip() or value.strip().lower() in {"nan", "none", "null", "na", "n/a"}
    try:
        return bool(pd.isna(value))
    except (TypeError, ValueError):
        return False


def _clean_text(value, fallback: str = "Unavailable") -> str:
    return fallback if _is_missing(value) else str(value).strip()


def _row_value(row: pd.Series, *names: str):
    for name in names:
        if name in row.index and not _is_missing(row[name]):
            return row[name]
    return None


def _count_from_gene_list(value) -> str:
    if _is_missing(value):
        return "Unavailable"
    if isinstance(value, (list, tuple, set)):
        return str(len(value))
    text = str(value).strip()
    try:
        parsed = json.loads(text)
        if isinstance(parsed, list):
            return str(len(parsed))
    except (json.JSONDecodeError, TypeError):
        pass
    return str(len([item for item in text.split(";") if item.strip()]))


def _format_count(value, gene_list_value=None) -> str:
    if _is_missing(value):
        return _count_from_gene_list(gene_list_value)
    text = str(value).strip()
    try:
        number = float(text)
        return str(int(number)) if math.isfinite(number) and number.is_integer() else text
    except (TypeError, ValueError):
        return text


def _format_fraction(value) -> str:
    if _is_missing(value):
        return "Unavailable"
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value).strip()
    return f"{number:.1%}" if math.isfinite(number) else "Unavailable"


def parse_warning_value(value) -> tuple[str, ...]:
    """Normalize list, JSON-list, semicolon, empty, and NaN warning values."""
    if _is_missing(value):
        return ()
    if isinstance(value, (list, tuple, set)):
        return tuple(str(item).strip() for item in value if not _is_missing(item))
    text = str(value).strip()
    try:
        parsed = json.loads(text)
        if isinstance(parsed, list):
            return tuple(str(item).strip() for item in parsed if not _is_missing(item))
        if isinstance(parsed, str):
            text = parsed.strip()
    except (json.JSONDecodeError, TypeError):
        try:
            parsed = ast.literal_eval(text)
            if isinstance(parsed, (list, tuple, set)):
                return tuple(str(item).strip() for item in parsed if not _is_missing(item))
        except (ValueError, SyntaxError):
            pass
    return tuple(part.strip() for part in text.split(";") if part.strip())


def classify_figure_filename(filename: str) -> tuple[str, str, int]:
    """Return human label, matched suffix, and display priority for a PNG name."""
    lower_name = filename.lower()
    for suffix, label in FIGURE_SUFFIX_LABELS:
        if lower_name.endswith(suffix.lower()):
            return label, suffix, FIGURE_PRIORITY[suffix]
    return Path(filename).stem, "", 100


def match_sample_from_figure(filename: str, samples: Iterable[str]) -> str | None:
    """Match a manifest sample prefix, preferring the longest sample name."""
    for sample in sorted({str(item) for item in samples if str(item)}, key=lambda item: (-len(item), item.lower())):
        if filename.lower().startswith(f"{sample}_".lower()):
            return sample
    return None


def sample_name_from_known_figure(filename: str) -> str | None:
    """Derive a sample only by removing a complete, known figure suffix."""
    lower_name = filename.lower()
    for suffix, _label in FIGURE_SUFFIX_LABELS:
        if lower_name.endswith(suffix.lower()):
            sample = filename[:-len(suffix)]
            return sample or None
    return None


def discover_figures(run_dir: Path, samples: Iterable[str]) -> dict[str, list[FigureRecord]]:
    sample_list = list(dict.fromkeys(str(sample) for sample in samples if str(sample)))
    found: dict[str, list[FigureRecord]] = {sample: [] for sample in sample_list}
    figures_dir = Path(run_dir) / "figures"
    if not figures_dir.is_dir():
        return found
    for path in figures_dir.iterdir():
        if not path.is_file() or path.suffix.lower() != ".png":
            continue
        sample = match_sample_from_figure(path.name, sample_list)
        if sample is None:
            continue
        label, suffix, priority = classify_figure_filename(path.name)
        if not suffix:
            sample_prefix = f"{sample}_"
            fallback = path.stem[len(sample_prefix):] if path.stem.lower().startswith(sample_prefix.lower()) else path.stem
            label = fallback or path.stem
        found[sample].append(FigureRecord(sample, path, label, suffix, priority))
    for sample in found:
        found[sample].sort(key=lambda item: (item.priority, item.label.lower(), item.path.name.lower()))
    return found


def _manifest_sample(row: pd.Series, row_number: int) -> str:
    sample = row.get("sample") if "sample" in row.index else None
    if not _is_missing(sample):
        return str(sample).strip()
    source = row.get("source_h5ad") if "source_h5ad" in row.index else None
    if not _is_missing(source):
        return Path(str(source)).stem
    return f"row_{row_number}"


def parse_manifest(manifest: pd.DataFrame | None) -> ManifestSummary:
    if manifest is None or not isinstance(manifest, pd.DataFrame):
        return ManifestSummary((), (), 0, "Run manifest unavailable or malformed.")
    if manifest.empty:
        return ManifestSummary((), (), 0, "Run manifest is empty.")
    successes: list[str] = []
    failures: list[tuple[str, str]] = []
    status_missing = "status" not in manifest.columns
    for position, (_, row) in enumerate(manifest.iterrows(), 1):
        sample = _manifest_sample(row, position)
        raw_status = row.get("status") if not status_missing else None
        status = _clean_text(raw_status, "status missing")
        if status.lower() == "ok":
            if sample not in successes:
                successes.append(sample)
        else:
            failures.append((sample, status))
    warning = "Run manifest has no status column." if status_missing else ""
    return ManifestSummary(tuple(successes), tuple(failures), len(manifest), warning)


def parse_context_summary_frame(table: pd.DataFrame) -> ContextSummary:
    if not isinstance(table, pd.DataFrame) or table.empty:
        return ContextSummary(False, (), "Context summary unavailable.")
    records: list[ContextFieldRecord] = []
    for _, row in table.iterrows():
        raw_field = _clean_text(_row_value(row, "field", "context_field"), "")
        normalized = raw_field.strip().upper().replace("_EXPR", "")
        if normalized not in {"H", "V"}:
            continue
        field = f"{normalized}_expr"
        requested_genes = _row_value(row, "requested_genes", "genes_requested")
        matched_genes = _row_value(row, "matched_genes", "genes_matched")
        warning_values: list[str] = []
        for column in ("warnings", "coverage_warning", "smoothing_warning"):
            if column in row.index:
                for warning in parse_warning_value(row[column]):
                    if warning not in warning_values:
                        warning_values.append(warning)
        records.append(ContextFieldRecord(
            field=field,
            status=_clean_text(_row_value(row, "status"), "Unavailable"),
            matched_gene_count=_format_count(
                _row_value(row, "matched_gene_count", "matched_count"), matched_genes
            ),
            requested_gene_count=_format_count(
                _row_value(row, "requested_gene_count", "requested_count"), requested_genes
            ),
            coverage_fraction=_format_fraction(
                _row_value(row, "coverage_fraction", "gene_coverage", "coverage")
            ),
            score_method=_clean_text(_row_value(row, "score_method", "method")),
            smoothing_method=_clean_text(_row_value(row, "smoothing_method", "smoothing")),
            high_state_fraction=_format_fraction(_row_value(row, "high_state_fraction")),
            warnings=tuple(warning_values),
        ))
    if not records:
        return ContextSummary(False, (), "Context summary unavailable: H_expr/V_expr rows were not found.")
    order = {"H_expr": 0, "V_expr": 1}
    records.sort(key=lambda item: order[item.field])
    return ContextSummary(True, tuple(records))


def load_context_summary(path: Path) -> ContextSummary:
    path = Path(path)
    if not path.is_file():
        return ContextSummary(False, (), "Context summary unavailable.")
    try:
        table = pd.read_csv(path, dtype=object, keep_default_na=True)
    except (OSError, UnicodeError, ValueError, pd.errors.ParserError, pd.errors.EmptyDataError) as exc:
        return ContextSummary(False, (), f"Context summary unavailable: {exc}")
    return parse_context_summary_frame(table)


def load_image_for_display(path: Path, width: int, height: int) -> tuple[Image.Image, tuple[int, int]]:
    """Read a PNG and return an aspect-preserving display copy and its original size."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError("The selected image file no longer exists.")
    with Image.open(path) as opened:
        opened.load()
        original_size = opened.size
        display = opened.convert("RGBA")
    display.thumbnail((max(1, int(width)), max(1, int(height))), Image.Resampling.LANCZOS)
    return display, original_size


def result_file_paths(run_dir: Path, sample: str) -> tuple[tuple[str, Path], ...]:
    run_dir = Path(run_dir)
    return (
        ("Manifest", run_dir / "spatial_graph_neighborhood_manifest.csv"),
        ("Cohort summary", run_dir / "combined_cohort_summary.csv"),
        ("Run parameters", run_dir / "run_parameters.json"),
        ("Context summary", run_dir / "context_fields" / f"{sample}_context_field_summary.csv"),
        ("Gene coverage", run_dir / "context_fields" / f"{sample}_context_field_gene_coverage.csv"),
        ("Leave-one-gene-out", run_dir / "context_fields" / f"{sample}_context_field_leave_one_gene_out.csv"),
        ("Neighborhood enrichment", run_dir / "neighborhood" / f"{sample}_categorical_enrichment.csv"),
        ("Binary association", run_dir / "neighborhood" / f"{sample}_binary_mask_association.csv"),
        ("Same-spot overlap", run_dir / "neighborhood" / f"{sample}_same_spot_overlap.csv"),
        ("Neighboring association", run_dir / "neighborhood" / f"{sample}_neighboring_spot_association.csv"),
        ("Continuous statistics", run_dir / "neighborhood" / f"{sample}_continuous_edge_statistics.csv"),
        ("Graph QC", run_dir / "spatial_graph" / f"{sample}_graph_qc.csv"),
        ("Graph metadata", run_dir / "spatial_graph" / f"{sample}_graph_metadata.json"),
        ("Degree distribution", run_dir / "spatial_graph" / f"{sample}_graph_degree_distribution.csv"),
    )


class SpatialGraphResultsPanel(ttk.Frame):
    """Read-only browser for the existing Spatial Graph output contract."""

    def __init__(self, parent, open_path: Callable[[Path], None]) -> None:
        super().__init__(parent)
        self.open_path = open_path
        self.run_dir: Path | None = None
        self.manifest_summary = ManifestSummary((), (), 0)
        self.figures_by_sample: dict[str, list[FigureRecord]] = {}
        self.current_figures: list[FigureRecord] = []
        self.current_figure_index = -1
        self.current_photo = None
        self.current_original_size: tuple[int, int] | None = None
        self._resize_job: str | None = None
        self.file_buttons: dict[str, ttk.Button] = {}
        self.sample_var = tk.StringVar()
        self.figure_var = tk.StringVar()
        self.status_var = tk.StringVar(value=INITIAL_MESSAGE)
        self.figure_info_var = tk.StringVar(value="No Spatial Graph run loaded.")
        self._build()
        self.refresh_control_states()

    def _build(self) -> None:
        ttk.Label(self, text="Spatial Graph Results", font=("Segoe UI Semibold", 13)).pack(anchor="w")
        ttk.Label(
            self,
            text="Read-only viewer for PNG, context QC, and result tables already written by Spatial Graph & Neighborhood.",
            foreground="#4b5563",
            wraplength=900,
        ).pack(anchor="w", pady=(2, 7))

        controls = ttk.Frame(self)
        controls.pack(fill="x")
        ttk.Label(controls, text="Sample").grid(row=0, column=0, sticky="w")
        ttk.Label(controls, text="Figure").grid(row=0, column=1, sticky="w", padx=(7, 0))
        self.sample_combo = ttk.Combobox(controls, textvariable=self.sample_var, state="readonly", width=24)
        self.sample_combo.grid(row=1, column=0, sticky="ew")
        self.figure_combo = ttk.Combobox(controls, textvariable=self.figure_var, state="readonly", width=48)
        self.figure_combo.grid(row=1, column=1, sticky="ew", padx=(7, 0))
        self.previous_button = ttk.Button(controls, text="Previous figure", command=lambda: self._step_figure(-1))
        self.previous_button.grid(row=1, column=2, padx=(7, 0))
        self.next_button = ttk.Button(controls, text="Next figure", command=lambda: self._step_figure(1))
        self.next_button.grid(row=1, column=3, padx=(5, 0))
        self.fit_button = ttk.Button(controls, text="Fit to window", command=self._display_current_figure)
        self.fit_button.grid(row=1, column=4, padx=(5, 0))
        self.open_image_button = ttk.Button(controls, text="Open image", command=self._open_current_image)
        self.open_image_button.grid(row=1, column=5, padx=(5, 0))
        self.open_folder_button = ttk.Button(controls, text="Open results folder", command=self._open_results_folder)
        self.open_folder_button.grid(row=1, column=6, padx=(5, 0))
        controls.columnconfigure(0, weight=1)
        controls.columnconfigure(1, weight=2)
        self.sample_combo.bind("<<ComboboxSelected>>", self._sample_changed)
        self.figure_combo.bind("<<ComboboxSelected>>", self._figure_changed)

        self.status_label = ttk.Label(self, textvariable=self.status_var, foreground="#4b5563", wraplength=930)
        self.status_label.pack(anchor="w", fill="x", pady=(6, 0))

        self.image_canvas = tk.Canvas(self, background="#111827", highlightthickness=0, height=330)
        self.image_canvas.pack(fill="both", expand=True, pady=(7, 0))
        self.image_canvas.bind("<Configure>", self._schedule_resize)
        self._show_placeholder(INITIAL_MESSAGE)
        ttk.Label(self, textvariable=self.figure_info_var, foreground="#4b5563", wraplength=930, justify="left").pack(
            anchor="w", fill="x", pady=(5, 0)
        )

        lower = ttk.Panedwindow(self, orient="horizontal")
        lower.pack(fill="x", pady=(7, 0))
        context_box = ttk.LabelFrame(lower, text="H_expr / V_expr context QC", padding=7)
        other_box = ttk.Frame(lower)
        lower.add(context_box, weight=3)
        lower.add(other_box, weight=2)
        self.context_text = tk.Text(context_box, height=10, wrap="word", state="disabled", background="#f8fafc")
        self.context_text.pack(fill="both", expand=True)
        self._set_text(self.context_text, "Context summary unavailable.")

        failures = ttk.LabelFrame(other_box, text="Failed or skipped samples", padding=7)
        failures.pack(fill="both", expand=True)
        self.failure_text = tk.Text(failures, height=4, wrap="word", state="disabled", background="#f8fafc")
        self.failure_text.pack(fill="both", expand=True)
        self._set_text(self.failure_text, "None recorded.")

        files = ttk.LabelFrame(other_box, text="Generated result files", padding=7)
        files.pack(fill="x", pady=(6, 0))
        file_labels = [label for label, _ in result_file_paths(Path("."), "sample")]
        for index, label in enumerate(file_labels):
            button = ttk.Button(files, text=label, command=lambda value=label: self._open_result_file(value))
            button.grid(row=index // 3, column=index % 3, sticky="ew", padx=(0 if index % 3 == 0 else 4, 0), pady=(0 if index < 3 else 4, 0))
            self.file_buttons[label] = button
        for column in range(3):
            files.columnconfigure(column, weight=1)

    @staticmethod
    def _set_text(widget: tk.Text, value: str) -> None:
        widget.configure(state="normal")
        widget.delete("1.0", "end")
        widget.insert("1.0", value)
        widget.configure(state="disabled")

    def load_run(self, run_dir: Path, manifest: pd.DataFrame) -> str:
        """Replace viewer state with one completed Spatial Graph run."""
        self.run_dir = Path(run_dir)
        self.current_photo = None
        self.current_original_size = None
        self.current_figures = []
        self.current_figure_index = -1
        self.sample_var.set("")
        self.figure_var.set("")
        self.sample_combo.configure(values=())
        self.figure_combo.configure(values=())
        self.manifest_summary = parse_manifest(manifest)
        samples = list(self.manifest_summary.successful_samples)
        placeholder_positions = [index for index, sample in enumerate(samples) if sample.startswith("row_")]
        if placeholder_positions:
            figures_dir = self.run_dir / "figures"
            inferred: list[str] = []
            if figures_dir.is_dir():
                for path in sorted(figures_dir.iterdir(), key=lambda item: item.name.lower()):
                    if path.is_file() and path.suffix.lower() == ".png":
                        sample = sample_name_from_known_figure(path.name)
                        if sample and sample not in inferred and sample not in samples:
                            inferred.append(sample)
            for index, sample in zip(placeholder_positions, inferred):
                samples[index] = sample
            self.manifest_summary = ManifestSummary(
                tuple(samples),
                self.manifest_summary.failed_samples,
                self.manifest_summary.total_rows,
                self.manifest_summary.warning,
            )
        self.figures_by_sample = discover_figures(self.run_dir, samples)
        self.sample_combo.configure(values=samples)
        self._display_failures()
        if not samples:
            message = "No successful samples were found in the run manifest."
            if self.manifest_summary.warning:
                message = f"{message} {self.manifest_summary.warning}"
            self.status_var.set(message)
            self.figure_info_var.set(f"Run directory: {self.run_dir}")
            self._show_placeholder(message)
            self._set_text(self.context_text, "Context summary unavailable.")
            self.refresh_control_states()
            return f"Spatial Graph & Neighborhood complete: 0/{self.manifest_summary.total_rows} samples. No successful Spatial Graph samples were available."
        self.sample_var.set(samples[0])
        self._load_selected_sample()
        figure_count = len(self.current_figures)
        if figure_count:
            self.status_var.set(f"Loaded {figure_count} figure(s) for sample {samples[0]}.")
            detail = f"Loaded {figure_count} figures for sample {samples[0]}."
        else:
            self.status_var.set("No Spatial Graph figures were generated. Review the run manifest.")
            detail = "No Spatial Graph figures were generated. Review the run manifest."
        return (
            f"Spatial Graph & Neighborhood complete: {len(samples)}/{self.manifest_summary.total_rows} samples. "
            f"{detail}"
        )

    def _display_failures(self) -> None:
        lines = [f"- {sample}: {status}" for sample, status in self.manifest_summary.failed_samples]
        if self.manifest_summary.warning:
            lines.append(f"- Manifest warning: {self.manifest_summary.warning}")
        self._set_text(self.failure_text, "\n".join(lines) if lines else "None recorded.")

    def _sample_changed(self, _event=None) -> None:
        self._load_selected_sample()

    def _load_selected_sample(self) -> None:
        sample = self.sample_var.get().strip()
        self.current_figures = list(self.figures_by_sample.get(sample, []))
        values = [item.selector_label for item in self.current_figures]
        self.figure_combo.configure(values=values)
        if self.current_figures:
            self.current_figure_index = 0
            self.figure_var.set(values[0])
            self._display_current_figure()
            self.status_var.set(f"Loaded {len(values)} figure(s) for sample {sample}.")
        else:
            self.current_figure_index = -1
            self.figure_var.set("")
            self._show_placeholder("No figures available for this sample.")
            self.figure_info_var.set(f"Sample: {sample}\nRun directory: {self.run_dir}\nNo figures available.")
            self.status_var.set("No figures available for this sample.")
        self._display_context_summary(sample)
        self.refresh_control_states()

    def _figure_changed(self, _event=None) -> None:
        selected = self.figure_var.get()
        for index, record in enumerate(self.current_figures):
            if record.selector_label == selected:
                self.current_figure_index = index
                break
        self._display_current_figure()
        self.refresh_control_states()

    def _step_figure(self, direction: int) -> None:
        if not self.current_figures:
            return
        self.current_figure_index = (self.current_figure_index + int(direction)) % len(self.current_figures)
        self.figure_var.set(self.current_figures[self.current_figure_index].selector_label)
        self._display_current_figure()
        self.refresh_control_states()

    def _schedule_resize(self, _event=None) -> None:
        if self._resize_job is not None:
            try:
                self.after_cancel(self._resize_job)
            except tk.TclError:
                pass
        self._resize_job = self.after(140, self._display_current_figure)

    def _show_placeholder(self, message: str) -> None:
        self.current_photo = None
        self.image_canvas.delete("all")
        width = max(self.image_canvas.winfo_width(), 320)
        height = max(self.image_canvas.winfo_height(), 180)
        self.image_canvas.create_text(
            width // 2,
            height // 2,
            text=message,
            fill="#cbd5e1",
            width=max(280, width - 60),
            justify="center",
        )

    def _display_current_figure(self) -> None:
        self._resize_job = None
        if not self.current_figures or not 0 <= self.current_figure_index < len(self.current_figures):
            return
        record = self.current_figures[self.current_figure_index]
        width = max(self.image_canvas.winfo_width() - 24, 100)
        height = max(self.image_canvas.winfo_height() - 24, 100)
        try:
            display, original_size = load_image_for_display(record.path, width, height)
            photo = ImageTk.PhotoImage(display, master=self.image_canvas)
        except FileNotFoundError:
            self.current_original_size = None
            self._show_placeholder("The selected image file no longer exists.")
            self.status_var.set("The selected image file no longer exists.")
            self.figure_info_var.set(f"Sample: {record.sample}\nFile: {record.path.name}\nRun directory: {self.run_dir}")
            self.refresh_control_states()
            return
        except (OSError, ValueError) as exc:
            self.current_original_size = None
            self._show_placeholder("The selected image could not be loaded.")
            self.status_var.set(f"The selected image could not be loaded: {exc}")
            self.figure_info_var.set(f"Sample: {record.sample}\nFile: {record.path.name}\nRun directory: {self.run_dir}")
            self.refresh_control_states()
            return
        self.current_photo = photo
        self.current_original_size = original_size
        self.image_canvas.delete("all")
        self.image_canvas.create_image(
            max(self.image_canvas.winfo_width(), 1) // 2,
            max(self.image_canvas.winfo_height(), 1) // 2,
            image=photo,
            anchor="center",
        )
        self.figure_info_var.set(
            f"Sample: {record.sample}\n"
            f"Figure: {record.label}\n"
            f"File: {record.path.name}\n"
            f"Size: {original_size[0]} × {original_size[1]}\n"
            f"Run directory: {self.run_dir}\n"
            f"Figure {self.current_figure_index + 1} of {len(self.current_figures)}"
        )

    def _display_context_summary(self, sample: str) -> None:
        if self.run_dir is None:
            self._set_text(self.context_text, "Context summary unavailable.")
            return
        path = self.run_dir / "context_fields" / f"{sample}_context_field_summary.csv"
        summary = load_context_summary(path)
        if not summary.available:
            self._set_text(self.context_text, summary.message or "Context summary unavailable.")
            return
        by_field = {record.field: record for record in summary.records}
        blocks: list[str] = []
        for field in ("H_expr", "V_expr"):
            record = by_field.get(field)
            if record is None:
                blocks.append(f"{field}\nStatus: Unavailable")
                continue
            status = "OK" if record.status.lower() == "ok" else (
                "Skipped by QC" if record.status.lower() == "skipped_qc" else record.status
            )
            warnings = "; ".join(record.warnings) if record.warnings else "None"
            blocks.append(
                f"{field}\n"
                f"Status: {status}\n"
                f"Genes: {record.matched_gene_count} / {record.requested_gene_count}\n"
                f"Coverage: {record.coverage_fraction}\n"
                f"Method: {record.score_method}\n"
                f"Smoothing: {record.smoothing_method}\n"
                f"High-state fraction: {record.high_state_fraction}\n"
                f"Warnings: {warnings}"
            )
        self._set_text(self.context_text, "\n\n".join(blocks))

    def _current_record(self) -> FigureRecord | None:
        if self.current_figures and 0 <= self.current_figure_index < len(self.current_figures):
            return self.current_figures[self.current_figure_index]
        return None

    def _safe_open(self, path: Path, label: str) -> None:
        if not path.exists():
            self.status_var.set(f"{label} is no longer available: {path}")
            self.refresh_control_states()
            return
        try:
            self.open_path(path)
        except Exception as exc:
            self.status_var.set(f"Could not open {label}: {exc}")
            messagebox.showerror(f"Open {label}", str(exc), parent=self)

    def _open_current_image(self) -> None:
        record = self._current_record()
        if record is not None:
            self._safe_open(record.path, "image")

    def _open_results_folder(self) -> None:
        if self.run_dir is not None:
            self._safe_open(self.run_dir, "results folder")

    def _available_result_files(self) -> dict[str, Path]:
        sample = self.sample_var.get().strip()
        if self.run_dir is None:
            return {}
        return {label: path for label, path in result_file_paths(self.run_dir, sample) if path.is_file()}

    def _open_result_file(self, label: str) -> None:
        path = self._available_result_files().get(label)
        if path is None:
            self.status_var.set(f"{label} is not available for the current run/sample.")
            self.refresh_control_states()
            return
        self._safe_open(path, label)

    def refresh_control_states(self) -> None:
        samples = tuple(self.sample_combo.cget("values"))
        figure_count = len(self.current_figures)
        current = self._current_record()
        self.sample_combo.configure(state="readonly" if samples else "disabled")
        self.figure_combo.configure(state="readonly" if figure_count else "disabled")
        navigation_state = "normal" if figure_count > 1 else "disabled"
        self.previous_button.configure(state=navigation_state)
        self.next_button.configure(state=navigation_state)
        self.fit_button.configure(state="normal" if current and current.path.is_file() else "disabled")
        self.open_image_button.configure(state="normal" if current and current.path.is_file() else "disabled")
        self.open_folder_button.configure(state="normal" if self.run_dir and self.run_dir.is_dir() else "disabled")
        available = self._available_result_files()
        for label, button in self.file_buttons.items():
            button.configure(state="normal" if label in available else "disabled")
