# SpatialTX Studio Desktop v0.2-beta

Released: 2026-07-05

v0.2-beta extends the publicly released v0.1-beta application. It does not replace or redefine the Transition Mapper, Cx, Sx, R(x), G(x), interface thresholds, or existing output files.

## Advanced Analysis

### Gene Composition

Reports every requested Cx and Sx gene, including missing genes, with mean and median transformed expression, detection fraction, and within-program relative contribution. Contribution percentages sum to 100% within each program when at least one present gene has nonzero weight.

### Interface Enrichment

Compares the unchanged v0.1 interface-like mask with all non-interface observations. Outputs include group means, composition percentages, fold enrichment, log2 fold enrichment, Hedges' g, two-sided Mann-Whitney p-values, and Benjamini-Hochberg FDR. Inferential statistics are left unavailable when group sizes are insufficient.

### Cx/Sx Interaction

Uses a six-nearest-neighbor graph and locally averaged, robustly scaled Cx/Sx fields. It reports coexistence, antagonism, balance, weighted spatial overlap, and Cx/Sx-dominant edge mixing. A seeded permutation null redistributes Sx activation over fixed coordinates and recomputes neighborhoods. Interaction results are not correlations.

## Output contract

Each Advanced Analysis run uses a new timestamped directory and creates a run manifest. Successful samples include CSV tables, 300-dpi PNG figures, vector PDF figures, and an `analysis_metadata.json` file recording definitions and parameters.

## In-app Results Dashboard

The **Run All 3 Analyses + Show Dashboard** action runs the three modules for every selected sample. The dashboard summarizes all module/sample results in one table. Selecting a summary row loads the complete corresponding CSV into the detailed table below, with direct actions for opening its CSV and publication figure. Individual module runs also update the same dashboard automatically.

The software remains a research prototype and is not intended for clinical decision-making.
