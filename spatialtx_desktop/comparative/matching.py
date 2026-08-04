from __future__ import annotations

import pandas as pd

from .models import SampleRecord


def pairwise_match(records: list[SampleRecord], reference: str, target: str) -> pd.DataFrame:
    reference_rows = [record for record in records if record.group == reference]
    target_rows = [record for record in records if record.group == target]
    if len(reference_rows) != 1 or len(target_rows) != 1 or len(records) != 2:
        raise ValueError("Pairwise comparison requires exactly one reference and one target sample.")
    return pd.DataFrame([{
        "comparison_id": f"{reference_rows[0].sample_id}_vs_{target_rows[0].sample_id}",
        "pair_id": "",
        "reference_sample_id": reference_rows[0].sample_id,
        "target_sample_id": target_rows[0].sample_id,
    }])


def paired_group_matches(records: list[SampleRecord], reference: str, target: str) -> pd.DataFrame:
    selected = [record for record in records if record.group in {reference, target}]
    missing_pair = [record.sample_id for record in selected if not record.pair_id]
    if missing_pair:
        raise ValueError(f"Paired mode requires pair_id for every selected sample: {', '.join(missing_pair)}")
    rows: list[dict] = []
    by_pair: dict[str, list[SampleRecord]] = {}
    for record in selected:
        by_pair.setdefault(record.pair_id, []).append(record)
    errors: list[str] = []
    for pair_id, members in sorted(by_pair.items()):
        ref = [record for record in members if record.group == reference]
        tar = [record for record in members if record.group == target]
        if len(ref) != 1 or len(tar) != 1 or len(members) != 2:
            errors.append(pair_id)
            continue
        rows.append({
            "comparison_id": pair_id,
            "pair_id": pair_id,
            "reference_sample_id": ref[0].sample_id,
            "target_sample_id": tar[0].sample_id,
        })
    if errors:
        raise ValueError(
            "Incomplete or duplicated reference/target members for pair_id value(s): " + ", ".join(errors)
        )
    if not rows:
        raise ValueError("No complete reference/target pairs were found.")
    return pd.DataFrame(rows)


def infer_manifest_batch_mode(records: list[SampleRecord], reference: str, target: str) -> str:
    selected = [record for record in records if record.group in {reference, target}]
    if selected and all(record.pair_id for record in selected):
        try:
            paired_group_matches(records, reference, target)
            return "paired"
        except ValueError:
            pass
    return "unpaired"
