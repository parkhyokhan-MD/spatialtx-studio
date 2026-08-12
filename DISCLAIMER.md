# Disclaimer

SpatialTX Studio Desktop v0.6-beta is an exploratory public source beta provided for research use only.

It is not a medical device and is not intended for diagnosis, prognosis, treatment selection, clinical decision-making, or any other clinical use. Outputs are computational summaries and hypothesis-generating candidates. They require independent review and validation before they can support scientific conclusions.

In particular, A3-A5 are optional exploratory utilities:

- A3 compares expression and detection summaries between two user-supplied conditions.
- A4 applies lightweight gene-symbol heuristics to prioritize receptor-like, membrane-associated, transporter-like, or surface-like candidates.
- A5 prepares a bounded candidate table for downstream QUBO-based combination selection.

These utilities do not discover or validate drug responses, receptor function, membrane localization, ligand-receptor binding, read-level evidence, biomarkers, biological subtypes, or clinical effects.

Multi-Pair Pre/Post results report observed spatial changes and a separate specimen comparability classification. `Good`, `Caution`, and `Low` describe transparent technical and sampling checks; they are not probabilities of biological truth. Regime transitions, directional arrows, and cohort counts must not be interpreted as confirmed treatment response, therapeutic efficacy, responder status, drug sensitivity, or predictive biomarker performance.

Users are responsible for input quality, study design, statistical analysis, interpretation, data governance, and compliance with applicable requirements. The software is provided without warranties or guarantees; see `LICENSE`.

## Development note

AI-assisted tools were used for documentation support, code organization, and troubleshooting during development. All scientific definitions, software behavior, release decisions, and final review were performed by the author.
