# SpatialTX Studio Desktop v0.6-dev

> Internal development source. Exploratory research use only; not for diagnosis, treatment selection, response prediction, or clinical decision-making.

## Positioning

v0.6-dev adds **Multiaxial Comparative Analysis + QC-aware paired interpretation** while preserving v0.5.5 and the primary C/S + FRAME2.6 model. It asks two separate questions: what changed between specimens, and how confidently the specimens can be compared.

## Added

- Optional paired H and V context axes using the existing SpatialTX programs.
- Visible Good/Caution/Low comparability beside regime and metric changes.
- QC-aware rule-based interpretation and explicit Low-comparability caution.
- `same_site`, `different_site`, and `unknown_site` metadata with SITE-SHIFT WARNING.
- Raw multiaxial Pre/Post/Delta export and a three-panel horizontal-delta figure.
- Compact primary mismatch summaries for available technical, sampling, and geometry QC.

New files:

- `context_changes.csv`
- `multiaxial_pair_summary.csv`
- `comparative_qc_summary.csv`
- `figures/multiaxial_pair_overview.png`

## Preserved boundaries

- `C(x)`, `S(x)`, and `R(x)=C(x)-S(x)` are unchanged.
- FRAME2.6 Type A/B/C definitions, Type B internal patterns, masks, thresholds, and transition metrics are unchanged.
- H/V do not modify C/S/R or operational regime assignment and may be unavailable without failing the pair.
- Existing Main Mapper, Single Pair, Import / Convert, GeoFlat, QUBO, Advanced Analysis, graph tools, and v0.5.5 exports remain available.
- No composite response score, direct treatment-effect attribution, deconvolution, new regime class, ML classifier, survival prediction, or comparative QUBO was added.

## Interpretation

H is a hypoxia-associated expression-context axis. V is an endothelial/angiogenic expression proxy, not perfusion, vessel density, or measured vascularity. All regime and directional outputs remain descriptive exploratory candidates requiring independent validation.
