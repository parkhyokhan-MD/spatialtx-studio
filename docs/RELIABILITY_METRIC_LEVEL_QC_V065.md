# SpatialTX Studio v0.65 Direction/CA_fraction metric-level QC

Status: internal development correction. Research use only; not for clinical decisions.

## Root-cause audit

The preserved score-domain implementation separates signed legacy Balance from nonnegative Activity. For the Activity path:

1. `valid_input` is true only when `activity_C` and `activity_S` are both finite and nonnegative.
2. `A = activity_C + activity_S` on `valid_input` spots.
3. `direction_defined = valid_input & (A > epsilon)`.
4. `Direction = (activity_C - activity_S) / (A + epsilon)` on the Direction-defined mask.
5. `CA_strength = A - abs(activity_C - activity_S)` on `valid_input` spots.
6. `CA_fraction = CA_strength / (A + epsilon)` where its explicit `ca_defined` mask is true.
7. In schema v3, `ca_defined` is stored independently even though it currently equals the Direction mask: `valid_input & (A > epsilon)`.

The old pair input gate counted `valid_input` spots, but the bootstrap and label-permutation helpers only required two finite values per group. Consequently, a globally valid pair could produce finite Direction or CA_fraction p-values from two defined spots. Those p-values entered `bh_adjust` in `build_pair_summary`. The BH helper itself was already correct: it adjusts finite values only and preserves `NaN` positions. The defect was the unsupported finite p-value generated before BH.

## Metric-level gate

Pair input validity remains separate. Direction and CA_fraction each export their own support fields:

- `defined_n`: `valid_input & metric_defined_mask & finite(metric_value)`.
- `valid_input_n`: number of Activity-valid spots.
- `defined_fraction = defined_n / valid_input_n`.
- QC status, inferential eligibility, and stable machine-readable reason.

Current masks and denominators:

```text
Direction defined_n:
sum(valid_input & direction_defined & finite(direction_D))

CA_fraction defined_n:
sum(valid_input & ca_defined & finite(ca_fraction))

denominator for both fractions:
valid_input_n
```

No all-spot denominator is used for the new metric-level gate. The historical all-spot Direction fraction remains available as `direction_defined_fraction_all_spots_legacy` in the internal summary for traceability.

Thresholds reuse the v0.65 validity settings:

| Status | Count | Defined fraction | Inference |
|---|---:|---:|---|
| PASS | `>= 30` | `>= 0.80` | eligible |
| CAUTION | `>= 30` | `>= 0.50` and `< 0.80` | not eligible |
| FAIL | `< 30`, or fraction `< 0.50`, or undefined/nonfinite metric | not eligible |

Both Pre and Post must be PASS for the pair-level Direction or CA_fraction inferential comparison. A sample with upstream Activity failure is FAIL; an upstream Activity warning prevents PASS. The pair itself is not invalidated solely because a derived metric lacks support.

## Inference and BH-FDR behavior

Descriptive medians and deltas remain available whenever mathematically defined. If Direction or CA_fraction is not inferentially eligible:

- bootstrap CI is `NaN`,
- permutation p-value is `NaN`,
- BH-FDR is `NaN`.

No placeholder `1.0` is generated. `bh_adjust` receives the complete metric family but operates only on finite eligible p-values; excluded positions remain `NaN`. Direction-based axis-dependence permutation rows use the same conservative sample support gate. Correlations may remain descriptive, but unsupported dependence p/FDR values are missing.

The preserved `delta_B_bh_fdr` compatibility field remains unchanged. The metric-QC-gated family is identified as Reliability schema v3; `delta_B_reliability_v2_bh_fdr` remains a compatibility alias and `delta_B_reliability_v3_bh_fdr` is the explicit current field.

## Exports and UI

Internal schema/cache identifier:

```text
v0.65-reliability-v3-metric-qc
```

New additive sidecar:

- `reliability_metric_qc.csv` — one row per pair, sample role, and axis.

The same QC status, eligibility, counts, fractions, reasons, denominator, and rule metadata are included where appropriate in `reliability_pair_summary.csv`, score-domain diagnostic CSV/JSON, `reliability_qc.json`, and run metadata. The six-pair validation tool also writes `six_pair_direction_ca_metric_qc.csv` and `.json`.

The v0.65 Reliability UI keeps descriptive values visible, shows Pre/Post defined support, and renders unsupported inference as `p: N/A — <reason>` rather than presenting `NaN` as a statistical result.

## Scope preserved

This correction does not change C, S, R, legacy Balance, interface/diffuse/transition burden, Type A/B/C, H/V, the Activity formula, the co-activation formula, or public v0.6 outputs. It adds no clipping, absolute-value conversion of inputs, shift, offset, min-max scaling, cohort correction, composite score, treatment-response claim, or clinical inference.
