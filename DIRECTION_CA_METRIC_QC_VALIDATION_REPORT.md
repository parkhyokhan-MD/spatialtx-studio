# SpatialTX Studio v0.65 Direction/CA_fraction metric-QC validation report

Validation date: 2026-08-24  
Status: PASS  
Scope: targeted metric-level support gate only; research use, non-diagnostic

## Confirmed root cause

The existing pair gate counted Activity `valid_input` spots, whereas Direction and CA_fraction are defined only where `A > epsilon`. The bootstrap and label-permutation helpers accepted two finite values per group. Therefore an otherwise valid pair could produce a finite Direction/CA_fraction p-value from only two defined spots. `bh_adjust` was already finite-only and preserved missing positions; the unsupported finite p-value was generated upstream in `build_pair_summary`.

## Implemented contract

```text
valid_input = finite(activity_C, activity_S)
              & activity_C >= 0
              & activity_S >= 0

A = activity_C + activity_S

direction_defined = valid_input & A > epsilon
D = (activity_C - activity_S) / (A + epsilon)

ca_defined = valid_input & A > epsilon
CA_strength = A - abs(activity_C - activity_S)
CA_fraction = CA_strength / (A + epsilon)

metric defined_n = sum(valid_input & metric_defined & finite(metric))
metric valid_input_n = sum(valid_input)
metric defined_fraction = defined_n / valid_input_n
```

Direction and CA currently share the same mathematical mask, but separate masks and explicit QC fields are retained for auditability.

| QC | Rule | Inferential eligibility |
|---|---|---|
| PASS | `defined_n >= 30` and `defined_fraction >= 0.80` | true |
| CAUTION | `defined_n >= 30` and `0.50 <= fraction < 0.80` | false |
| FAIL | `defined_n < 30`, fraction `< 0.50`, no defined spots, or nonfinite metric | false |

Both Pre and Post must be PASS. If not eligible, descriptive values remain available where finite, while bootstrap CI, permutation p, and BH-FDR are `NaN`/null. Missing p-values are excluded from BH and remain missing afterward.

## Automated tests

- Full suite: **203 passed, 1 skipped, 26 subtests passed** in 138.34 seconds.
- New cases cover 2 spots, 29 spots, the inclusive 30-count and 80%-fraction boundaries, 79%, exactly 50%, 49%, zero defined spots, finite-only BH exclusion, explicit export fields, Direction dependence gating, and GUI `p: N/A` rendering.
- Exact `defined_n=30` and `defined_fraction=0.80` cannot occur simultaneously because `30 / 0.80 = 37.5` spots. The two inclusive boundaries are therefore tested separately: `30/36` for the count boundary and `32/40` for the exact 80% boundary.

## v0.6 Balance regression

- Result: **72/72 PASS**.
- Tolerance: `rtol=atol=1e-12`.
- Preserved fields: Pre/Post/delta Balance, Balance CI/p/compatibility FDR, and C/S gene coverage.
- Existing C/S/R and v0.6-compatible outputs did not change.

## Six-pair real-data metric QC

All 12 slides passed Activity input validity. `defined_fraction` uses `valid_input_n` as denominator.

| Pair | Valid input Pre/Post | Direction defined Pre/Post | Direction fraction Pre/Post | Direction QC | Eligible | CA_fraction defined Pre/Post | CA fraction Pre/Post | CA QC | Eligible |
|---|---:|---:|---:|---|---|---:|---:|---|---|
| sample_30 | 431 / 145 | 370 / 71 | 85.85% / 48.97% | FAIL | false | 370 / 71 | 85.85% / 48.97% | FAIL | false |
| sample_38 | 453 / 245 | 352 / 123 | 77.70% / 50.20% | CAUTION | false | 352 / 123 | 77.70% / 50.20% | CAUTION | false |
| sample_41 | 536 / 1040 | 448 / 868 | 83.58% / 83.46% | PASS | true | 448 / 868 | 83.58% / 83.46% | PASS | true |
| sample_42 | 143 / 996 | 138 / 969 | 96.50% / 97.29% | PASS | true | 138 / 969 | 96.50% / 97.29% | PASS | true |
| sample_43 | 694 / 886 | 340 / 828 | 48.99% / 93.45% | FAIL | false | 340 / 828 | 48.99% / 93.45% | FAIL | false |
| sample_44 | 1274 / 548 | 969 / 507 | 76.06% / 92.52% | CAUTION | false | 969 / 507 | 76.06% / 92.52% | CAUTION | false |

Pair counts for both Direction and CA_fraction: **PASS 2, CAUTION 2, FAIL 2**. Ineligible four pairs retained missing Direction/CA_fraction p/FDR values; only samples 41 and 42 ran those inferential comparisons.

## Seed-42 determinism

Two independent full six-pair runs produced identical SHA-256 values for **18/18** Reliability artifacts, including:

- spot, pair-summary, metric-QC, coverage, exclusivity, dependence, and score-domain CSV files;
- Reliability diagnostic/QC and six-pair metric-QC JSON files;
- dependence and valid-fraction PNG files;
- regression fixed summary/diff/report artifacts.

`run_metadata.json` was excluded because it intentionally contains the run tag and creation timestamp. No deterministic computation field was excluded.

## Remaining limits

- The resampling unit is the spot distribution of unregistered slides, not a biological replicate; no treatment-effect or specimen-level claim is supported.
- The Activity source is log1p count-like abundance without library-size normalization and can remain sensitive to depth, composition, and occupancy.
- The 30/80/50 thresholds are predefined QC rules with boundary discontinuities; they are not learned or biologically calibrated.
- Descriptive Direction/CA_fraction values with low support can be unstable even though they are retained for audit.
- Non-targeted Balance, Activity, and CA_strength inference rules remain unchanged and should be reviewed separately before any stronger inferential interpretation.
- H/V remain observational context axes, and there is no composite response score.

## Internal schema and new outputs

Internal schema: `v0.65-reliability-v3-metric-qc`.

New primary output: `reliability_metric_qc.csv`. The pair summary, score-domain diagnostic CSV/JSON, Reliability QC JSON, and run metadata also expose the relevant counts, fractions, statuses, eligibility, rules, and reasons. The six-pair tool adds `six_pair_direction_ca_metric_qc.csv` and `.json`.

No GitHub commit, push, release, or remote modification was performed.

## Source audit

Preserved baseline ZIP SHA-256:

```text
100A2AC09F1524D0C13779BDBF99FAC3113F125FDE56302FC50F0A2E9A7B8782
```

Changed files compared with that ZIP:

- `CHANGELOG.md`
- `README.md`
- `README_DESKTOP.md`
- `docs/RELIABILITY_LAYER_V065.md`
- `docs/RELIABILITY_SCORE_DOMAIN_FIX_V065.md`
- `spatialtx_desktop/comparative/multi_pair.py`
- `spatialtx_desktop/multi_pair_ui.py`
- `spatialtx_desktop/reliability/core.py`
- `spatialtx_desktop/reliability/dependence.py`
- `spatialtx_desktop/reliability/diagnostics.py`
- `spatialtx_desktop/reliability/exports.py`
- `spatialtx_desktop/reliability/models.py`
- `tests/test_comparative_ui.py`
- `tests/test_reliability_layer.py`
- `tools/validate_v06_six_pair.py`

New files:

- `DIRECTION_CA_METRIC_QC_VALIDATION_REPORT.md`
- `docs/RELIABILITY_METRIC_LEVEL_QC_V065.md`

No baseline file was removed.
