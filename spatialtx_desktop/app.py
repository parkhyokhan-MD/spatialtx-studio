from __future__ import annotations

import os
import queue
import threading
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, ttk

from . import APP_NAME, AUTHOR, BUILD_DATE, __version__
from .advanced_ui import AdvancedToolsPanel
from .workflow import (
    DEFAULT_C_GENES,
    DEFAULT_S_GENES,
    export_result,
    optimize_genes,
    parse_gene_text,
    run_batch,
    save_optimizer_results,
    scan_h5ad,
)


class SpatialTXDesktop(tk.Tk):
    def __init__(self) -> None:
        super().__init__()
        self.title(f"{APP_NAME} v{__version__}")
        width = min(1680, max(1180, int(self.winfo_screenwidth() * .90)))
        height = min(940, max(760, int(self.winfo_screenheight() * .86)))
        self.geometry(f"{width}x{height}")
        self.minsize(1180, 720)
        self.option_add("*Font", ("Segoe UI", 9))
        self.files: list[Path] = []
        self.last_run: Path | None = None
        self.last_summary = None
        self.pending_genes: dict[str, list[str] | None] = {"C": None, "S": None}
        self.advanced_qubo_pool: list[str] | None = None
        self.map_photo = None
        self.map_paths: list[Path] = []
        self.busy = False
        self.messages: queue.Queue[tuple[str, object]] = queue.Queue()
        self._build_style()
        self._build_ui()
        self.after(100, self._poll_messages)

    def _build_style(self) -> None:
        style = ttk.Style(self)
        if "vista" in style.theme_names():
            style.theme_use("vista")
        style.configure("Title.TLabel", font=("Segoe UI Semibold", 18))
        style.configure("Sub.TLabel", foreground="#4b5563")
        style.configure("Primary.TButton", font=("Segoe UI Semibold", 10))
        style.configure("Treeview", rowheight=24)

    def _build_ui(self) -> None:
        header = ttk.Frame(self, padding=(18, 14, 18, 8))
        header.pack(fill="x")
        ttk.Label(header, text="SpatialTX Studio", style="Title.TLabel").pack(side="left")
        ttk.Label(header, text=f"Desktop v{__version__}  •  h5ad workflow  •  {AUTHOR}", style="Sub.TLabel").pack(side="left", padx=12, pady=(7, 0))
        ttk.Label(header, text="Research use only — not for clinical decisions", foreground="#9a3412").pack(side="right", pady=(7, 0))

        main = ttk.Panedwindow(self, orient="horizontal")
        main.pack(fill="both", expand=True, padx=16, pady=(0, 12))
        left = ttk.Frame(main, padding=2)
        right = ttk.Frame(main, padding=(12, 2, 2, 2))
        main.add(left, weight=2)
        main.add(right, weight=3)

        self.right_tabs = ttk.Notebook(right)
        self.right_tabs.pack(fill="both", expand=True)
        analysis_tab = ttk.Frame(self.right_tabs, padding=10)
        map_tab = ttk.Frame(self.right_tabs, padding=10)
        qubo_tab = ttk.Frame(self.right_tabs, padding=10)
        theory_tab = ttk.Frame(self.right_tabs, padding=10)
        interpretation_tab = ttk.Frame(self.right_tabs, padding=10)
        about_tab = ttk.Frame(self.right_tabs, padding=10)
        advanced_tab = ttk.Frame(self.right_tabs, padding=10)
        self.right_tabs.add(analysis_tab, text="Analysis")
        self.right_tabs.add(map_tab, text="Map Viewer")
        self.right_tabs.add(qubo_tab, text="QUBO Optimizer")
        self.right_tabs.add(theory_tab, text="Theory & Metrics")
        self.right_tabs.add(interpretation_tab, text="Interpretation")
        self.right_tabs.add(advanced_tab, text="Advanced Tools")
        self.right_tabs.add(about_tab, text="About & Version")
        self.interpretation_tab = interpretation_tab

        source = ttk.LabelFrame(left, text="1  Input and samples", padding=10)
        source.pack(fill="both", expand=True)
        folder_row = ttk.Frame(source)
        folder_row.pack(fill="x")
        self.input_var = tk.StringVar()
        ttk.Entry(folder_row, textvariable=self.input_var).pack(side="left", fill="x", expand=True)
        ttk.Button(folder_row, text="Browse…", command=self._browse_input).pack(side="left", padx=(6, 0))
        ttk.Button(folder_row, text="Scan h5ad", command=self._scan).pack(side="left", padx=(6, 0))

        self.sample_tree = ttk.Treeview(source, columns=("sample", "size", "path"), show="headings", selectmode="extended", height=11)
        self.sample_tree.heading("sample", text="Sample")
        self.sample_tree.heading("size", text="Size")
        self.sample_tree.heading("path", text="Path")
        self.sample_tree.column("sample", width=180, anchor="w")
        self.sample_tree.column("size", width=72, anchor="e")
        self.sample_tree.column("path", width=420, anchor="w")
        sample_scroll = ttk.Scrollbar(source, orient="vertical", command=self.sample_tree.yview)
        self.sample_tree.configure(yscrollcommand=sample_scroll.set)
        self.sample_tree.pack(side="left", fill="both", expand=True, pady=(8, 0))
        sample_scroll.pack(side="right", fill="y", pady=(8, 0))

        sample_actions = ttk.Frame(left, padding=(10, 4, 10, 0))
        sample_actions.pack(fill="x")
        ttk.Button(sample_actions, text="Select all", command=lambda: self.sample_tree.selection_set(self.sample_tree.get_children())).pack(side="left")
        ttk.Button(sample_actions, text="Clear", command=lambda: self.sample_tree.selection_remove(self.sample_tree.selection())).pack(side="left", padx=5)
        self.sample_count = ttk.Label(sample_actions, text="No files scanned", style="Sub.TLabel")
        self.sample_count.pack(side="right")

        programs = ttk.LabelFrame(left, text="2  C/S gene programs", padding=10)
        programs.pack(fill="x", pady=(10, 0))
        gene_grid = ttk.Frame(programs)
        gene_grid.pack(fill="x")
        ttk.Label(gene_grid, text="C-side (immune)").grid(row=0, column=0, sticky="w")
        ttk.Label(gene_grid, text="S-side (stromal)").grid(row=0, column=1, sticky="w", padx=(10, 0))
        self.c_text = tk.Text(gene_grid, height=4, wrap="word", undo=True)
        self.s_text = tk.Text(gene_grid, height=4, wrap="word", undo=True)
        self.c_text.grid(row=1, column=0, sticky="nsew")
        self.s_text.grid(row=1, column=1, sticky="nsew", padx=(10, 0))
        gene_grid.columnconfigure(0, weight=1); gene_grid.columnconfigure(1, weight=1)
        self.c_text.insert("1.0", ", ".join(DEFAULT_C_GENES))
        self.s_text.insert("1.0", ", ".join(DEFAULT_S_GENES))

        settings = ttk.LabelFrame(analysis_tab, text="3  Scoring and output", padding=10)
        settings.pack(fill="x")
        qrow = ttk.Frame(settings); qrow.pack(fill="x")
        self.cq, self.sq, self.gq = tk.StringVar(value="0.80"), tk.StringVar(value="0.80"), tk.StringVar(value="0.60")
        for index, (label, var) in enumerate((("C quantile", self.cq), ("S quantile", self.sq), ("G quantile", self.gq))):
            box = ttk.Frame(qrow); box.grid(row=0, column=index, sticky="ew", padx=(0 if index == 0 else 7, 0))
            ttk.Label(box, text=label).pack(anchor="w")
            ttk.Entry(box, textvariable=var, width=9).pack(fill="x")
            qrow.columnconfigure(index, weight=1)
        outrow = ttk.Frame(settings); outrow.pack(fill="x", pady=(10, 0))
        self.output_var = tk.StringVar(value=str(Path.cwd() / "desktop_results"))
        ttk.Entry(outrow, textvariable=self.output_var).pack(side="left", fill="x", expand=True)
        ttk.Button(outrow, text="Output…", command=self._browse_output).pack(side="left", padx=(6, 0))
        ttk.Button(settings, text="Run SpatialTX scoring + maps", style="Primary.TButton", command=self._run).pack(fill="x", pady=(10, 0))
        ttk.Label(
            settings,
            text="Quantiles define relative high-score spots within each sample. See Theory & Metrics before changing defaults.",
            style="Sub.TLabel", wraplength=400,
        ).pack(fill="x", pady=(6, 0))

        results = ttk.LabelFrame(analysis_tab, text="4  Results and export", padding=10)
        results.pack(fill="both", expand=True, pady=(10, 0))
        result_buttons = ttk.Frame(results); result_buttons.pack(fill="x")
        ttk.Button(result_buttons, text="Open results", command=self._open_results).pack(side="left")
        ttk.Button(result_buttons, text="Export folder…", command=lambda: self._export(False)).pack(side="left", padx=5)
        ttk.Button(result_buttons, text="Export zip…", command=lambda: self._export(True)).pack(side="left")
        self.log = tk.Text(results, height=12, state="disabled", wrap="word", background="#f8fafc")
        self.log.pack(fill="both", expand=True, pady=(8, 0))
        self.progress = ttk.Progressbar(results, mode="indeterminate")
        self.progress.pack(fill="x", pady=(7, 0))

        self._build_theory_tab(theory_tab)
        self._build_interpretation_tab(interpretation_tab)
        self._build_about_tab(about_tab)
        self._build_qubo_tab(qubo_tab)
        self._build_map_tab(map_tab)
        AdvancedToolsPanel(advanced_tab, on_qubo_pool=self._set_advanced_qubo_pool).pack(fill="both", expand=True)

    def _build_qubo_tab(self, parent: ttk.Frame) -> None:
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
        ttk.Label(parent, text="C/S gene-program optimizer", font=("Segoe UI Semibold", 13)).pack(anchor="w")
        ttk.Label(
            parent,
            text="Select exactly one h5ad sample on the left. Optimize and apply each side independently, then recompute maps when ready.",
            style="Sub.TLabel", wraplength=440,
        ).pack(anchor="w", pady=(2, 10))

        settings = ttk.LabelFrame(parent, text="1  Optimization settings", padding=10)
        settings.pack(fill="x")
        optrow = ttk.Frame(settings); optrow.pack(fill="x")
        self.k_var, self.pool_var = tk.IntVar(value=8), tk.IntVar(value=40)
        self.iter_var, self.seed_var = tk.IntVar(value=300), tk.IntVar(value=20260624)
        for index, (label, var) in enumerate((("Genes (k)", self.k_var), ("Pool", self.pool_var), ("Iterations", self.iter_var), ("Seed", self.seed_var))):
            box = ttk.Frame(optrow); box.grid(row=0, column=index, sticky="ew", padx=(0 if index == 0 else 7, 0))
            ttk.Label(box, text=label).pack(anchor="w")
            upper = 99999999 if label == "Seed" else 2000 if label == "Iterations" else 100
            ttk.Spinbox(box, from_=1 if label == "Seed" else 2, to=upper, textvariable=var, width=9).pack(fill="x")
            optrow.columnconfigure(index, weight=1)
        guide = ttk.LabelFrame(settings, text="What do these options mean?", padding=8)
        guide.pack(fill="x", pady=(9, 0))
        ttk.Label(
            guide,
            text=(
                "Genes (k): the fixed number of genes selected for the optimized program. Smaller values are compact and easier to "
                "interpret; larger values retain broader signal but may add redundancy. Default: 8.\n\n"
                "Pool: the maximum number of candidate genes considered. A larger pool broadens the search but can increase runtime "
                "and variability. Default: 40.\n\n"
                "Iterations: simulated-annealing swap attempts. More iterations can improve the search at the cost of runtime. "
                "This does not change k. Default: 300.\n\n"
                "Seed: fixes the random search path so the same input and settings reproduce the same result. Default: 20260624."
            ),
            wraplength=410, justify="left",
        ).pack(anchor="w", fill="x")
        self.qubo_pool_status = ttk.Label(settings, text="Candidate source: built-in side markers + variable genes", style="Sub.TLabel", wraplength=410)
        self.qubo_pool_status.pack(anchor="w", fill="x", pady=(6, 0))
        method = ttk.LabelFrame(settings, text="How QUBO optimization works", padding=8)
        method.pack(fill="x", pady=(8, 0))
        ttk.Label(
            method,
            text=(
                "1. Build a bounded candidate gene pool.\n"
                "2. Score each gene for C/S alignment, directional R signal, gradient association, spatial enrichment, detection, and variance.\n"
                "3. Penalize overlap with the opposite side, low detection, and highly redundant gene pairs.\n"
                "4. Convert the task into a binary optimization problem that must select exactly k genes.\n"
                "5. Solve it quickly on the local CPU with a classical simulated-annealing heuristic.\n"
                "6. Apply the selected program, then recompute C/S fields and redraw the maps.\n\n"
                "Key idea: QUBO does not simply rank genes one at a time. It selects a complementary gene combination that explains "
                "the C-side or S-side spatial direction while avoiding unnecessary redundancy."
            ),
            wraplength=620, justify="left",
        ).pack(anchor="w", fill="x")

        sides = ttk.Frame(parent)
        sides.pack(fill="both", expand=True, pady=(10, 0))
        sides.columnconfigure(0, weight=1); sides.columnconfigure(1, weight=1)
        for column, side in enumerate(("C", "S")):
            title = "2  C-side (immune)" if side == "C" else "3  S-side (stromal)"
            frame = ttk.LabelFrame(sides, text=title, padding=10)
            frame.grid(row=0, column=column, sticky="nsew", padx=(0, 5) if side == "C" else (5, 0))
            ttk.Button(frame, text=f"Run QUBO for {side}-side", command=lambda value=side: self._optimize(value)).pack(fill="x")
            ttk.Label(frame, text="Selected genes", style="Sub.TLabel").pack(anchor="w", pady=(8, 2))
            selected_text = tk.Text(frame, height=8, wrap="word", state="disabled", background="#f8fafc")
            selected_text.pack(fill="both", expand=True)
            apply_button = ttk.Button(frame, text=f"Apply {side} genes only", command=lambda value=side: self._apply_optimized_side(value), state="disabled")
            apply_button.pack(fill="x", pady=(7, 0))
            ttk.Button(frame, text=f"Restore fixed {side} gene set", command=lambda value=side: self._restore_fixed_side(value)).pack(fill="x", pady=(5, 0))
            status = ttk.Label(frame, text="Not optimized yet.", style="Sub.TLabel", wraplength=190)
            status.pack(fill="x", pady=(5, 0))
            setattr(self, f"qubo_{side.lower()}_text", selected_text)
            setattr(self, f"apply_{side.lower()}_button", apply_button)
            setattr(self, f"qubo_{side.lower()}_status", status)

        finish = ttk.LabelFrame(parent, text="4  Recompute and redraw", padding=10)
        finish.pack(fill="x", pady=(10, 0))
        ttk.Label(
            finish,
            text="Apply C, S, or both above. This final action runs the currently displayed C/S programs on all samples selected at left.",
            style="Sub.TLabel", wraplength=420,
        ).pack(anchor="w")
        ttk.Button(finish, text="Restore both fixed gene sets", command=self._restore_both_fixed).pack(fill="x", pady=(7, 0))
        ttk.Button(finish, text="Recompute scoring + redraw maps", style="Primary.TButton", command=self._recompute_from_qubo).pack(fill="x", pady=(7, 0))

    def _build_about_tab(self, parent: ttk.Frame) -> None:
        ttk.Label(parent, text=APP_NAME, font=("Segoe UI Semibold", 18), foreground="#0f3d56").pack(anchor="w", pady=(8, 0))
        ttk.Label(parent, text=f"Version {__version__}", font=("Segoe UI Semibold", 12)).pack(anchor="w", pady=(2, 14))

        card = ttk.LabelFrame(parent, text="Build information", padding=14)
        card.pack(fill="x")
        details = [
            ("Creator", AUTHOR),
            ("Version", f"v{__version__}"),
            ("Release date", BUILD_DATE),
            ("Edition", "First public beta release for Windows"),
        ]
        for row, (label, value) in enumerate(details):
            ttk.Label(card, text=label, font=("Segoe UI Semibold", 9)).grid(row=row, column=0, sticky="nw", pady=4)
            ttk.Label(card, text=value, wraplength=310).grid(row=row, column=1, sticky="nw", padx=(14, 0), pady=4)
        card.columnconfigure(1, weight=1)

        description = ttk.LabelFrame(parent, text="Version description", padding=14)
        description.pack(fill="both", expand=True, pady=(12, 0))
        ttk.Label(
            description,
            text=(
                f"SpatialTX Studio Desktop v{__version__} is the first public Windows research prototype. "
                "It scores immune-side C(x) and stromal-side S(x) gene programs, constructs the balance field "
                "R(x)=C(x)-S(x), estimates its local spatial gradient G(x), and generates exploratory interface and diffuse-transition "
                "calls. The desktop edition also provides a classical fixed-cardinality QUBO-inspired optimizer for compact C/S gene programs.\n\n"
                "Operational regime labels are exploratory candidates and are not validated biological subtypes."
            ),
            wraplength=430, justify="left",
        ).pack(anchor="nw", fill="x")
        ttk.Label(
            description,
            text="Research use only — not intended for diagnosis, treatment selection, or clinical decision-making.",
            foreground="#9a3412", wraplength=430, justify="left",
        ).pack(anchor="w", pady=(16, 0))

    def _build_map_tab(self, parent: ttk.Frame) -> None:
        toolbar = ttk.Frame(parent)
        toolbar.pack(fill="x")
        ttk.Label(toolbar, text="Spatial map").pack(side="left")
        self.map_sample_var = tk.StringVar()
        self.map_combo = ttk.Combobox(toolbar, textvariable=self.map_sample_var, state="readonly", width=28)
        self.map_combo.pack(side="left", padx=6)
        self.map_combo.bind("<<ComboboxSelected>>", lambda _event: self._display_selected_map())
        ttk.Button(toolbar, text="Previous", command=lambda: self._step_map(-1)).pack(side="left")
        ttk.Button(toolbar, text="Next", command=lambda: self._step_map(1)).pack(side="left", padx=5)
        ttk.Button(toolbar, text="Fit to window", command=self._display_selected_map).pack(side="left")
        ttk.Button(toolbar, text="Open externally", command=self._open_viewer_map).pack(side="right")
        self.map_canvas = tk.Canvas(parent, background="#111827", highlightthickness=0)
        self.map_canvas.pack(fill="both", expand=True, pady=(8, 0))
        self.map_canvas.create_text(20, 20, anchor="nw", fill="#cbd5e1", text="Run SpatialTX scoring to view generated maps here.", tags="placeholder")
        self.map_canvas.bind("<Configure>", lambda _event: self.after_idle(self._display_selected_map) if self.map_paths else None)

    def _refresh_map_list(self) -> None:
        self.map_paths = []
        names: list[str] = []
        if self.last_summary is not None:
            for _, row in self.last_summary.iterrows():
                path = Path(str(row.get("spatial_map_png", "")))
                if str(row.get("status", "")) == "ok" and path.is_file():
                    self.map_paths.append(path)
                    names.append(str(row.get("sample", path.stem)))
        self.map_combo.configure(values=names)
        if names:
            self.map_combo.current(0)
            self._display_selected_map()

    def _display_selected_map(self) -> None:
        if not self.map_paths or not hasattr(self, "map_canvas"):
            return
        from PIL import Image, ImageTk

        index = self.map_combo.current()
        if index < 0 or index >= len(self.map_paths):
            index = 0
        try:
            with Image.open(self.map_paths[index]) as source:
                image = source.convert("RGB")
                width = max(200, self.map_canvas.winfo_width() - 20)
                height = max(200, self.map_canvas.winfo_height() - 20)
                image.thumbnail((width, height), Image.Resampling.LANCZOS)
                self.map_photo = ImageTk.PhotoImage(image)
            self.map_canvas.delete("all")
            self.map_canvas.create_image(self.map_canvas.winfo_width() // 2, self.map_canvas.winfo_height() // 2, image=self.map_photo, anchor="center")
        except Exception as exc:
            self.map_canvas.delete("all")
            self.map_canvas.create_text(20, 20, anchor="nw", fill="#fecaca", text=f"Could not display map:\n{exc}")

    def _step_map(self, direction: int) -> None:
        if not self.map_paths:
            return
        self.map_combo.current((max(0, self.map_combo.current()) + direction) % len(self.map_paths))
        self._display_selected_map()

    def _open_viewer_map(self) -> None:
        index = self.map_combo.current()
        if not self.map_paths or index < 0:
            messagebox.showinfo("Spatial map", "Run scoring first.", parent=self); return
        try:
            os.startfile(self.map_paths[index])
        except OSError as exc:
            messagebox.showerror("Spatial map", str(exc), parent=self)

    def _set_advanced_qubo_pool(self, genes: list[str], source: Path) -> None:
        self.advanced_qubo_pool = list(genes)
        self.pool_var.set(max(len(genes), self.k_var.get()))
        self.qubo_pool_status.configure(text=f"Candidate source: Advanced export ({len(genes)} genes) — {source}")
        self._write_log(f"Loaded {len(genes)} Advanced candidates into the QUBO optimizer.")

    def _build_theory_tab(self, parent: ttk.Frame) -> None:
        toolbar = ttk.Frame(parent)
        toolbar.pack(fill="x")
        ttk.Label(toolbar, text="Interface-centered model", font=("Segoe UI Semibold", 13)).pack(side="left")
        ttk.Button(toolbar, text="Back to analysis", command=lambda: self.right_tabs.select(0)).pack(side="right")

        container = ttk.Frame(parent)
        container.pack(fill="both", expand=True, pady=(8, 0))
        scroll = ttk.Scrollbar(container, orient="vertical")
        theory = tk.Text(
            container, wrap="word", yscrollcommand=scroll.set, relief="flat",
            background="#f8fafc", padx=14, pady=12, spacing1=3, spacing3=7,
        )
        scroll.configure(command=theory.yview)
        theory.pack(side="left", fill="both", expand=True)
        scroll.pack(side="right", fill="y")
        theory.tag_configure("h1", font=("Segoe UI Semibold", 14), foreground="#0f3d56", spacing1=10, spacing3=5)
        theory.tag_configure("h2", font=("Segoe UI Semibold", 11), foreground="#155e75", spacing1=8, spacing3=3)
        theory.tag_configure("formula", font=("Consolas", 10), background="#eaf4f7", lmargin1=12, lmargin2=12, spacing1=5, spacing3=5)
        theory.tag_configure("note", foreground="#7c2d12", background="#fff7ed", lmargin1=8, lmargin2=8)
        sections = [
            ("h1", "What SpatialTX is asking\n"),
            (None, "SpatialTX does not classify a tissue by cell type alone. It asks whether immune-side and stromal-side expression programs form a localized spatial interface, a diffuse transition field, or a relatively flat pattern.\n"),
            ("h2", "1. Core fields\n"),
            ("formula", "C(x) = mean z-scored expression of the C-side gene program\nS(x) = mean z-scored expression of the S-side gene program\nR(x) = C(x) - S(x)\nG(x) = mean |R(x) - R(neighbor)| over six nearest spatial neighbors\n"),
            (None, "C and S are relative program activities within a sample. R is the local balance: positive values favor the C-side program and negative values favor the S-side program. G becomes large where that balance changes sharply in space.\n"),
            ("h2", "2. Interface and diffuse transition calls\n"),
            ("formula", "high_C = C >= sample C quantile\nhigh_S = S >= sample S quantile\nhigh_G = G >= sample G quantile\ninterface = high_C AND high_S AND high_G\ndiffuse = high_G AND (high_C OR high_S) AND NOT interface\n"),
            (None, "The defaults (C=.80, S=.80, G=.60) are sample-relative thresholds. Raising a quantile makes the call more selective; lowering it increases sensitivity and usually increases the number of candidate spots. These are exploratory settings, not universal biological cutoffs.\n"),
            ("h2", "3. Operational regimes\n"),
            (None, "Type_A_candidate: interface fraction is at least 1%; a localized interface-like signal is present.\n\nType_B_candidate: Type A is absent, but diffuse fraction is at least 5%; transition burden is spatially distributed rather than localized.\n\nType_C_candidate: neither condition is met; the sample is transition-poor or comparatively flat under the current gene programs and thresholds.\n"),
            ("h2", "4. Metrics worth reading together\n"),
            (None, "interface_fraction — fraction of spots meeting the interface rule.\ninterface_coherence_score — interface fraction multiplied by the largest connected-component ratio; favors a coherent interface over scattered hits.\ndiffuse_fraction — fraction of diffuse-transition spots.\ntransition_burden_score — combines diffuse fraction, upper-tail gradient strength, and spatial fragmentation.\nR_dynamic_range — separation between the 90th and 10th percentiles of R.\nC/S_gene_coverage — requested program genes actually present in the h5ad feature set. Low coverage weakens interpretation.\nQC_flag — PASS, WARN, or FAIL based on gene coverage, coordinate validity, feature-name uniqueness, C/S overlap, and spot count.\n"),
            ("h2", "5. QUBO-inspired optimizer\n"),
            (None, "The optimizer chooses exactly k genes from a bounded candidate pool. It rewards agreement with the selected side, directional R separation, spatial-gradient alignment, variance, and detection. It penalizes opposite-side correlation, low detection, and redundant gene pairs. A classical simulated-annealing search minimizes the resulting binary objective. It is QUBO-inspired, not a quantum backend.\n"),
            ("note", "Interpretation guardrail\nRegime labels are operational candidates. Gene-program edits and optimization can change the result, so compare baseline and recomputed maps and review gene coverage. SpatialTX Studio is an exploratory research prototype and is not intended for diagnosis or treatment decisions.\n"),
        ]
        for tag, text in sections:
            theory.insert("end", text, tag or ())
        theory.configure(state="disabled")

    def _build_interpretation_tab(self, parent: ttk.Frame) -> None:
        ttk.Label(parent, text="Run interpretation", font=("Segoe UI Semibold", 13)).pack(anchor="w")
        ttk.Label(
            parent,
            text="After scoring, this table links each operational call to the metrics that produced it.",
            style="Sub.TLabel", wraplength=430,
        ).pack(anchor="w", pady=(2, 8))
        columns = ("sample", "regime", "qc", "interface", "coherence", "diffuse", "burden")
        self.result_tree = ttk.Treeview(parent, columns=columns, show="headings", height=8, selectmode="browse")
        headings = {"sample": "Sample", "regime": "Regime", "qc": "QC", "interface": "Interface", "coherence": "Coherence", "diffuse": "Diffuse", "burden": "Burden"}
        widths = {"sample": 120, "regime": 125, "qc": 48, "interface": 68, "coherence": 72, "diffuse": 62, "burden": 62}
        for column in columns:
            self.result_tree.heading(column, text=headings[column])
            self.result_tree.column(column, width=widths[column], anchor="w" if column in {"sample", "regime", "qc"} else "e")
        self.result_tree.pack(fill="x")
        self.result_tree.bind("<<TreeviewSelect>>", self._on_result_selected)
        actions = ttk.Frame(parent)
        actions.pack(fill="x", pady=(6, 0))
        ttk.Button(actions, text="Open selected map", command=self._open_selected_map).pack(side="left")
        ttk.Button(actions, text="Open results folder", command=self._open_results).pack(side="left", padx=6)
        self.interpret_text = tk.Text(parent, wrap="word", height=16, state="disabled", background="#f8fafc", padx=10, pady=8)
        self.interpret_text.pack(fill="both", expand=True, pady=(8, 0))
        self._set_interpretation_text("No completed run yet. Run SpatialTX scoring, then return here for sample-level interpretation.")

    def _browse_input(self) -> None:
        folder = filedialog.askdirectory(title="Select folder containing h5ad files")
        if folder:
            self.input_var.set(folder)
            self._scan()

    def _browse_output(self) -> None:
        folder = filedialog.askdirectory(title="Select results folder")
        if folder:
            self.output_var.set(folder)

    def _scan(self) -> None:
        try:
            self.files = scan_h5ad(self.input_var.get())
        except Exception as exc:
            messagebox.showerror("Scan failed", str(exc), parent=self)
            return
        self.sample_tree.delete(*self.sample_tree.get_children())
        for index, path in enumerate(self.files):
            size = path.stat().st_size / (1024 * 1024)
            self.sample_tree.insert("", "end", iid=str(index), values=(path.stem, f"{size:.1f} MB", str(path)))
        self.sample_count.configure(text=f"{len(self.files)} h5ad file(s)")
        self._write_log(f"Scanned {len(self.files)} h5ad file(s) in {self.input_var.get()}")

    def _selected(self) -> list[Path]:
        return [self.files[int(item)] for item in self.sample_tree.selection()]

    def _genes(self) -> tuple[list[str], list[str]]:
        c = parse_gene_text(self.c_text.get("1.0", "end"))
        s = parse_gene_text(self.s_text.get("1.0", "end"))
        if not c or not s:
            raise ValueError("Both C-side and S-side gene programs must contain at least one gene.")
        return c, s

    def _quantiles(self) -> tuple[float, float, float]:
        values = tuple(float(v.get()) for v in (self.cq, self.sq, self.gq))
        if any(v <= 0 or v >= 1 for v in values):
            raise ValueError("C, S and G quantiles must be between 0 and 1.")
        return values

    def _background(self, function) -> None:
        self.busy = True
        self.progress.start(12)
        self._set_controls(False)
        threading.Thread(target=function, daemon=True).start()

    def _set_controls(self, enabled: bool) -> None:
        # Keep the UI readable while a data job is active; command handlers also validate state.
        self.config(cursor="" if enabled else "watch")

    def _run(self) -> None:
        if self.busy:
            messagebox.showinfo("SpatialTX is busy", "Wait for the current job to finish.", parent=self); return
        try:
            selected, (c, s), quantiles = self._selected(), self._genes(), self._quantiles()
            if not selected:
                raise ValueError("Select at least one sample in the scan table.")
            output = self.output_var.get().strip()
            if not output:
                raise ValueError("Choose an output folder.")
        except Exception as exc:
            messagebox.showwarning("Cannot run", str(exc), parent=self); return

        def work():
            try:
                run_dir, summary = run_batch(selected, output, c, s, lambda m: self.messages.put(("log", m)), *quantiles)
                self.messages.put(("run_done", (run_dir, summary)))
            except Exception as exc:
                self.messages.put(("error", ("SpatialTX run failed", exc)))
        self._background(work)

    def _optimize(self, side: str) -> None:
        if self.busy:
            messagebox.showinfo("SpatialTX is busy", "Wait for the current job to finish.", parent=self); return
        try:
            selected, (c, s) = self._selected(), self._genes()
            if len(selected) != 1:
                raise ValueError("Select exactly one sample for gene optimization.")
            k, pool, iterations, seed = self.k_var.get(), self.pool_var.get(), self.iter_var.get(), self.seed_var.get()
            if k > pool:
                raise ValueError("Selected gene count (k) cannot exceed candidate pool size.")
        except Exception as exc:
            messagebox.showwarning("Cannot optimize", str(exc), parent=self); return
        self._write_log(f"Optimizing {side}-side genes for {selected[0].name}…")

        def work():
            try:
                genes, detail, summary = optimize_genes(
                    selected[0], side, c, s, k, pool, iterations,
                    candidate_genes=self.advanced_qubo_pool, seed=seed,
                )
                root = self.last_run or Path(self.output_var.get())
                saved = save_optimizer_results(root, selected[0].stem, side, detail, summary)
                self.messages.put(("optimizer_done", (side, genes, summary, saved)))
            except Exception as exc:
                self.messages.put(("error", ("Optimizer failed", exc)))
        self._background(work)

    def _apply_optimized_side(self, side: str) -> None:
        genes = self.pending_genes.get(side)
        if not genes:
            messagebox.showinfo("QUBO optimizer", f"Run the {side}-side optimizer first.", parent=self); return
        target = self.c_text if side == "C" else self.s_text
        target.delete("1.0", "end"); target.insert("1.0", ", ".join(genes))
        status = self.qubo_c_status if side == "C" else self.qubo_s_status
        status.configure(text=f"Applied to the {side} gene box. Maps have not been recomputed yet.")
        self._write_log(f"Applied optimized {side}-side program only: {', '.join(genes)}")

    def _restore_fixed_side(self, side: str) -> None:
        genes = DEFAULT_C_GENES if side == "C" else DEFAULT_S_GENES
        target = self.c_text if side == "C" else self.s_text
        target.delete("1.0", "end")
        target.insert("1.0", ", ".join(genes))
        self.pending_genes[side] = None
        text_widget = self.qubo_c_text if side == "C" else self.qubo_s_text
        apply_button = self.apply_c_button if side == "C" else self.apply_s_button
        status = self.qubo_c_status if side == "C" else self.qubo_s_status
        text_widget.configure(state="normal")
        text_widget.delete("1.0", "end")
        text_widget.insert("1.0", "Fixed set restored:\n" + "\n".join(genes))
        text_widget.configure(state="disabled")
        apply_button.configure(state="disabled")
        status.configure(text=f"Fixed {side} program restored. Recompute maps when ready.")
        self._write_log(f"Restored fixed {side}-side gene program: {', '.join(genes)}")

    def _restore_both_fixed(self) -> None:
        self._restore_fixed_side("C")
        self._restore_fixed_side("S")
        self._write_log("Both fixed C/S gene programs are restored. Use recompute/redraw to generate fixed-set results.")

    def _recompute_from_qubo(self) -> None:
        self._write_log("Recomputing and redrawing with the currently applied C/S gene programs.")
        self._run()

    @staticmethod
    def _number(row, key: str, default: float = 0.0) -> float:
        try:
            value = float(row.get(key, default))
            return value if value == value else default
        except (TypeError, ValueError):
            return default

    def _render_interpretation(self, summary) -> None:
        self.result_tree.delete(*self.result_tree.get_children())
        for index, row in summary.reset_index(drop=True).iterrows():
            status = str(row.get("status", ""))
            if status != "ok":
                values = (row.get("sample", ""), "ERROR", "FAIL", "-", "-", "-", "-")
            else:
                values = (
                    row.get("sample", ""), row.get("regime_label", ""), row.get("QC_flag", ""),
                    f"{self._number(row, 'interface_fraction'):.3f}",
                    f"{self._number(row, 'interface_coherence_score'):.3f}",
                    f"{self._number(row, 'diffuse_fraction'):.3f}",
                    f"{self._number(row, 'transition_burden_score'):.3f}",
                )
            self.result_tree.insert("", "end", iid=str(index), values=values)
        children = self.result_tree.get_children()
        if children:
            self.result_tree.selection_set(children[0])
            self.result_tree.focus(children[0])
            self._show_interpretation_row(0)

    def _on_result_selected(self, _event=None) -> None:
        selection = self.result_tree.selection()
        if selection:
            self._show_interpretation_row(int(selection[0]))

    def _show_interpretation_row(self, index: int) -> None:
        if self.last_summary is None or index >= len(self.last_summary):
            return
        row = self.last_summary.reset_index(drop=True).iloc[index]
        status = str(row.get("status", ""))
        if status != "ok":
            self._set_interpretation_text(
                f"{row.get('sample', 'Sample')} could not be interpreted because scoring failed.\n\n{status}\n\n"
                "Review the run log, gene presence, and spatial coordinates in the input h5ad."
            )
            return
        regime = str(row.get("regime_label", "Unclassified"))
        regime_text = {
            "Type_A_candidate": (
                "A localized interface-like candidate was detected. High C, high S, and high local balance-gradient spots "
                "co-occur above the current sample-relative thresholds. Inspect whether these spots form a plausible tissue boundary."
            ),
            "Type_B_candidate": (
                "A localized interface was not prominent, but diffuse transition burden was detected. The balance field changes across "
                "multiple regions rather than concentrating into one interface."
            ),
            "Type_C_candidate": (
                "The sample is transition-poor or comparatively flat under the current C/S programs and thresholds. This does not prove "
                "biological absence; low feature coverage or unsuitable programs can produce the same operational call."
            ),
        }.get(regime, "No standard operational regime was assigned.")
        interface = self._number(row, "interface_fraction")
        coherence = self._number(row, "interface_coherence_score")
        diffuse = self._number(row, "diffuse_fraction")
        burden = self._number(row, "transition_burden_score")
        dynamic = self._number(row, "R_dynamic_range")
        c_cov = self._number(row, "C_gene_coverage")
        s_cov = self._number(row, "S_gene_coverage")
        coverage_note = "Gene-program coverage is adequate for a first-pass read."
        if min(c_cov, s_cov) < .5:
            coverage_note = "CAUTION: at least one program has <50% feature coverage. Treat the biological interpretation as weak."
        elif min(c_cov, s_cov) < .8:
            coverage_note = "Some requested genes are missing. Compare results after reviewing the detected feature names."
        pattern = str(row.get("public_transition_pattern", "") or "")
        pattern_line = f"\nDiffuse pattern: {pattern}" if pattern and pattern.lower() != "nan" else ""
        qc_flag = str(row.get("QC_flag", "PASS"))
        qc_notes = str(row.get("QC_notes", "") or "")
        qc_line = f"• QC: {qc_flag}" + (f" ({qc_notes})" if qc_notes and qc_notes.lower() != "nan" else "")
        text = (
            f"{row.get('sample', 'Sample')} — {regime}\n\n{regime_text}\n{pattern_line}\n\n"
            "Key measurements\n"
            f"• Interface fraction: {interface:.3f}\n"
            f"• Interface coherence: {coherence:.3f}\n"
            f"• Diffuse fraction: {diffuse:.3f}\n"
            f"• Transition burden: {burden:.3f}\n"
            f"• R dynamic range: {dynamic:.3f}\n"
            f"• C/S gene coverage: {c_cov:.0%} / {s_cov:.0%}\n\n"
            f"{qc_line}\n\n"
            f"{coverage_note}\n\n"
            "Suggested review\n"
            "1. Open the six-panel map and compare C(x), S(x), R(x), and G(x).\n"
            "2. Check whether candidate spots follow tissue morphology rather than isolated noise.\n"
            "3. Review missing program genes before changing thresholds.\n"
            "4. If using optimized genes, compare baseline and recomputed maps; do not interpret optimization as validation.\n\n"
            "Operational exploratory result — not a clinical classification."
        )
        self._set_interpretation_text(text)

    def _set_interpretation_text(self, text: str) -> None:
        self.interpret_text.configure(state="normal")
        self.interpret_text.delete("1.0", "end")
        self.interpret_text.insert("1.0", text)
        self.interpret_text.configure(state="disabled")

    def _open_selected_map(self) -> None:
        selection = self.result_tree.selection()
        if not selection or self.last_summary is None:
            messagebox.showinfo("Spatial map", "Select a completed sample first.", parent=self); return
        row = self.last_summary.reset_index(drop=True).iloc[int(selection[0])]
        path = Path(str(row.get("spatial_map_png", "")))
        if not path.is_file():
            messagebox.showerror("Spatial map", f"Map file was not found:\n{path}", parent=self); return
        try:
            os.startfile(path)
        except OSError as exc:
            messagebox.showerror("Spatial map", str(exc), parent=self)

    def _open_results(self) -> None:
        if not self.last_run:
            messagebox.showinfo("Results", "Run scoring first.", parent=self); return
        try:
            os.startfile(self.last_run)
        except OSError as exc:
            messagebox.showerror("Results", str(exc), parent=self)

    def _export(self, as_zip: bool) -> None:
        if not self.last_run:
            messagebox.showinfo("Export", "Run scoring first.", parent=self); return
        destination = (filedialog.asksaveasfilename(title="Export results zip", defaultextension=".zip", filetypes=[("ZIP archive", "*.zip")])
                       if as_zip else filedialog.askdirectory(title="Export results folder into…"))
        if not destination:
            return
        try:
            exported = export_result(self.last_run, destination, as_zip)
            self._write_log(f"Exported: {exported}")
            messagebox.showinfo("Export complete", str(exported), parent=self)
        except Exception as exc:
            messagebox.showerror("Export failed", str(exc), parent=self)

    def _poll_messages(self) -> None:
        try:
            while True:
                kind, payload = self.messages.get_nowait()
                if kind == "log":
                    self._write_log(str(payload))
                elif kind == "run_done":
                    self.last_run, summary = payload
                    self.last_summary = summary
                    self._render_interpretation(summary)
                    self._refresh_map_list()
                    ok = int((summary["status"] == "ok").sum()) if "status" in summary else 0
                    self._write_log(f"Run finished: {ok}/{len(summary)} successful\n{self.last_run}")
                    self.progress.stop(); self.busy = False; self._set_controls(True)
                    self.right_tabs.select(self.interpretation_tab)
                    messagebox.showinfo("SpatialTX complete", f"Results saved to:\n{self.last_run}", parent=self)
                elif kind == "optimizer_done":
                    side, genes, summary, saved = payload
                    self.pending_genes[side] = genes
                    text_widget = self.qubo_c_text if side == "C" else self.qubo_s_text
                    apply_button = self.apply_c_button if side == "C" else self.apply_s_button
                    status = self.qubo_c_status if side == "C" else self.qubo_s_status
                    text_widget.configure(state="normal")
                    text_widget.delete("1.0", "end")
                    text_widget.insert("1.0", "\n".join(genes))
                    text_widget.configure(state="disabled")
                    apply_button.configure(state="normal")
                    status.configure(text=f"{len(genes)} genes selected. Review, then apply this side.")
                    self._write_log(f"Optimizer selected {len(genes)} {side}-side genes. Details: {saved}")
                    self.progress.stop(); self.busy = False; self._set_controls(True)
                elif kind == "error":
                    title, exc = payload
                    self._write_log(f"ERROR: {exc}")
                    self.progress.stop(); self.busy = False; self._set_controls(True)
                    messagebox.showerror(title, str(exc), parent=self)
        except queue.Empty:
            pass
        self.after(100, self._poll_messages)

    def _write_log(self, message: str) -> None:
        self.log.configure(state="normal")
        self.log.insert("end", str(message).rstrip() + "\n")
        self.log.see("end")
        self.log.configure(state="disabled")


def main() -> None:
    app = SpatialTXDesktop()
    app.mainloop()


if __name__ == "__main__":
    main()
