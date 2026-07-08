# SpatialTX Studio Desktop developer notes

Release version: v0.3-beta

Release date: 2026-07-08

## Current scope

- Adds a dedicated **Import / Convert** tab.
- Converts supported Raw 10x MEX/MTX and Raw Visium H5 + spatial folders to validated AnnData `.h5ad` files.
- Keeps Advanced / Experimental focused on hypothesis-generation analysis rather than input-format conversion.
- Keeps `.h5ad` as the only canonical input analyzed by the Main Mapper.
- Adds lightweight robustness and memory-safety diagnostics without changing conservative defaults.
- Does not add RDS, h5Seurat, parquet, generic CSV, or new biological interpretation support.

## Development guardrails

- Preserve the established Main Mapper defaults unless a future release explicitly documents a breaking change.
- Keep raw-format conversion separate from core H5AD analysis.
- Keep Advanced / Experimental outputs clearly labeled as hypothesis-generation or review utilities.
- Do not represent operational Type A/B/C labels, receptor-like filtering, QUBO selection, or threshold stability as biological validation or clinical evidence.
