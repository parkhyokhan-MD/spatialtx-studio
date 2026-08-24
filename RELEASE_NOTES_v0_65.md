# SpatialTX Studio Desktop v0.65

Release date: 2026-08-24  
Research use only; not for clinical decisions.

## Highlights

- Adds the default-off Multi-axis Reliability Layer without changing established v0.6 C/S/R, spatial masks, H/V context, comparability, or Type A/B/C candidate outputs.
- Preserves signed legacy Balance and uses a separate nonnegative pre-z-score program abundance for Activity, Direction, and Co-activation.
- Adds Direction and CA_fraction metric-level support QC using `defined_n / valid_input_n`.
- Requires at least 30 defined spots and at least 80% defined support in both Pre and Post before Direction/CA_fraction CI, permutation p-value, or BH-FDR is calculated.
- Retains descriptive values when possible while representing unsupported inference as `NaN`/null in exports and `p: N/A — reason` in the GUI.
- Adds `reliability_metric_qc.csv`, stable QC reasons, score-domain diagnostics, finite-eligible-only BH handling, and internal schema `v0.65-reliability-v3-metric-qc`.

## Validation

- Final v0.65 packaging suite: 204 passed, including the previously environment-dependent Tk test; 26 Reliability subtests passed.
- Preserved Balance regression: 72/72 comparisons PASS at `rtol=atol=1e-12`.
- Six-pair Direction QC: PASS 2, CAUTION 2, FAIL 2.
- Six-pair CA_fraction QC: PASS 2, CAUTION 2, FAIL 2.
- Two independent seed-42 runs: 18/18 Reliability artifacts byte-identical by SHA-256, excluding only timestamp/run-tag metadata.

See `DIRECTION_CA_METRIC_QC_VALIDATION_REPORT.md`, `docs/RELIABILITY_LAYER_V065.md`, and `docs/RELIABILITY_METRIC_LEVEL_QC_V065.md` for formulas, audit details, edge-case tests, and limitations.

## Included v0.65 screenshots

- `docs/screenshots/spatialtx_studio_v0_65_reliability_metric_qc_overview.png`
- `docs/screenshots/spatialtx_studio_v0_65_score_domain_audit.png`
- `docs/screenshots/spatialtx_studio_v0_65_spatial_context_joint_state.png`

The H/V joint-state map is observational context only. No registration, spot-wise subtraction, composite response score, specimen-level treatment inference, or clinical prediction is performed.
