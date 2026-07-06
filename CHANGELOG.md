# Changelog

All notable public changes to SpatialTX Studio Desktop are recorded here.

## 0.2-beta — 2026-07-05

- Preserved the v0.1-beta Transition Mapper workflow, Cx/Sx definitions, defaults, and output contracts.
- Added an **Advanced Analysis** workspace with Gene Composition, Interface Enrichment, and Cx/Sx Interaction tabs.
- Added per-gene relative-contribution tables and manuscript-ready bar charts.
- Added interface/non-interface composition, fold enrichment, Hedges' g, Mann-Whitney testing, and BH-FDR.
- Added neighborhood-based coexistence, antagonism, balance, weighted spatial overlap, edge mixing, and seeded permutation inference.
- Added 300-dpi PNG, vector PDF, CSV, run-manifest, and JSON provenance outputs.
- Added a separate `advanced_cli.py` entry point; the original `app_cli.py` behavior is unchanged.
- Added an in-app Results Dashboard with a combined module/sample summary, full CSV table viewer, direct CSV/figure actions, and a one-click run-all workflow.
- Expanded Theory & Metrics with the rationale, formulas, statistical assumptions, spatial-permutation reference, and interpretation limitations for all three Advanced Analysis modules.

## 0.1-beta — 2026-07-02

- First public source release.
- Added local desktop workflows for `.h5ad` discovery, spatial scoring, map generation, interpretation, and export.
- Added local C-side and S-side gene-program optimization with a classical simulated-annealing fallback.
- Added opt-in exploratory utilities for MEX conversion, condition comparison, heuristic candidate filtering, and QUBO candidate-pool handoff.
- Added QC flags, gene-coverage reporting, coordinate fallbacks, and research-use guardrails.
- Added source installation guidance, release documentation, citation metadata, and licensing notices.
