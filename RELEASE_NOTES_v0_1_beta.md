# SpatialTX Studio Desktop v0.1-beta

Release date: 2026-07-02

This is the first public source beta of the local desktop research prototype. It supports exploratory analysis of spatial transcriptomics `.h5ad` files, C/S program scoring, transition-map generation, QC review, and result export.

## Important limitations

- Research use only; not for diagnosis, treatment selection, or clinical decision-making.
- Results are exploratory and depend on input preprocessing, spatial coordinates, gene coverage, thresholds, and selected programs.
- Operational regime labels are candidates, not validated biological subtypes.
- QUBO-based selection is a local, classical optimization aid and does not establish biological validity.
- A3-A5 are hypothesis-generation utilities. They do not discover or validate drug responses, receptors, ligand-receptor interactions, biomarkers, or clinical effects.
- Advanced raw-data and sequence utilities are local conversion or template helpers; they do not query or validate against external biological databases.

## Release package

`SpatialTX_Studio_Desktop_v0_1_beta_SOURCE.zip` contains the Python source, launch scripts, examples, configuration, and release documentation. It does not contain datasets, generated results, caches, build products, or local analysis outputs.

Review `README.md`, `README_DESKTOP.md`, `DISCLAIMER.md`, and `THIRD_PARTY_LICENSES.md` before use.
