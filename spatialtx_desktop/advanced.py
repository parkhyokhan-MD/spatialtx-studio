from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd


KNOWN_RECEPTORS = {
    "EGFR", "ERBB2", "ERBB3", "ERBB4", "MET", "AXL", "PDGFRA", "PDGFRB", "FGFR1", "FGFR2",
    "FGFR3", "KDR", "FLT1", "TGFBR1", "TGFBR2", "IL6R", "IL6ST", "CXCR3", "CXCR4", "CCR2",
    "CCR5", "PDCD1", "CD274", "CTLA4", "TIGIT", "HAVCR2", "LAG3", "FAS", "TNFRSF1A",
}
KNOWN_LIGANDS = {
    "EGF", "TGFA", "HGF", "AREG", "EREG", "VEGFA", "VEGFB", "PDGFA", "PDGFB", "FGF1", "FGF2",
    "TGFB1", "TGFB2", "IL6", "IL10", "IFNG", "TNF", "CXCL9", "CXCL10", "CXCL12", "CCL2", "CCL5",
    "CD274", "FASLG", "SPP1", "GAS6", "WNT5A",
}


def scan_pre_post_pairs(folder: str | Path) -> pd.DataFrame:
    files = sorted(Path(folder).rglob("*.h5ad"))
    groups: dict[str, dict[str, Path]] = {}
    for path in files:
        stem = path.stem
        match = re.search(r"(?i)(?:^|[_\-.])(pre|post)(?:[_\-.]|$)", stem)
        if not match:
            continue
        stage = match.group(1).lower()
        key = re.sub(r"(?i)(?:^|[_\-.])(pre|post)(?:[_\-.]|$)", "_", stem).strip("_.-")
        groups.setdefault(key or "sample", {})[stage] = path.resolve()
    return pd.DataFrame([
        {"pair": key, "pre_h5ad": str(value["pre"]), "post_h5ad": str(value["post"]), "status": "paired"}
        for key, value in groups.items() if "pre" in value and "post" in value
    ])


def _mean_detection(adata) -> tuple[np.ndarray, np.ndarray, str]:
    X = adata.X
    raw_values = np.asarray(X.data if hasattr(X, "data") and not isinstance(X, np.ndarray) else X).ravel()
    finite_values = raw_values[np.isfinite(raw_values)]
    count_like = bool(
        finite_values.size
        and np.nanmin(finite_values) >= 0
        and np.mean(np.isclose(finite_values, np.round(finite_values), atol=1e-6)) >= .98
    )
    if hasattr(X, "getnnz"):
        values = X.tocsr().astype(float)
        detection = np.asarray(X.getnnz(axis=0)).ravel() / max(1, adata.n_obs)
        if count_like:
            totals = np.asarray(values.sum(axis=1)).ravel()
            scale = np.divide(10000.0, totals, out=np.zeros_like(totals, dtype=float), where=totals > 0)
            values = values.multiply(scale[:, None]).tocsr()
            values.data = np.log1p(values.data)
        mean = np.asarray(values.mean(axis=0)).ravel()
    else:
        values = np.asarray(X, dtype=float)
        detection = np.mean(values > 0, axis=0)
        if count_like:
            totals = np.nansum(values, axis=1)
            scale = np.divide(10000.0, totals, out=np.zeros_like(totals, dtype=float), where=totals > 0)
            values = np.log1p(values * scale[:, None])
        mean = np.nanmean(values, axis=0)
    return mean, detection, "normalize_total_1e4_log1p" if count_like else "existing_processed_scale"


def compare_pre_post(pre_path: str | Path, post_path: str | Path, output_csv: str | Path) -> pd.DataFrame:
    import anndata as ad

    pre, post = ad.read_h5ad(pre_path), ad.read_h5ad(post_path)
    pre_mean, pre_detection, pre_scale = _mean_detection(pre)
    post_mean, post_detection, post_scale = _mean_detection(post)
    pre_lookup = {str(g).upper(): i for i, g in enumerate(pre.var_names)}
    post_lookup = {str(g).upper(): i for i, g in enumerate(post.var_names)}
    if len(pre_lookup) != pre.n_vars or len(post_lookup) != post.n_vars:
        raise ValueError("Pre/post comparison requires unique feature names, ignoring letter case.")
    common = sorted(pre_lookup.keys() & post_lookup.keys())
    if not common:
        raise ValueError("Pre and post h5ad files have no shared feature names.")
    rows = []
    for key in common:
        i, j = pre_lookup[key], post_lookup[key]
        delta = float(post_mean[j] - pre_mean[i])
        delta_detection = float(post_detection[j] - pre_detection[i])
        logfc_like = float(np.log2((post_mean[j] + 1e-6) / (pre_mean[i] + 1e-6)))
        direction = "post_up" if delta > 0 else "pre_up" if delta < 0 else "no_change"
        candidate_score = abs(delta) * (0.5 + max(float(pre_detection[i]), float(post_detection[j]))) + 0.25 * abs(delta_detection)
        rows.append({
            "gene": str(post.var_names[j]), "mean_pre": float(pre_mean[i]), "mean_post": float(post_mean[j]),
            "delta_mean": delta, "logfc_like": logfc_like,
            "detection_pre": float(pre_detection[i]), "detection_post": float(post_detection[j]),
            "delta_detection": delta_detection, "direction": direction, "candidate_score": candidate_score,
            "pre_expression_scale": pre_scale, "post_expression_scale": post_scale,
        })
    result = pd.DataFrame(rows).sort_values("candidate_score", ascending=False)
    output = Path(output_csv); output.parent.mkdir(parents=True, exist_ok=True); result.to_csv(output, index=False)
    return result


def filter_receptor_membrane(input_csv: str | Path, output_csv: str | Path) -> pd.DataFrame:
    table = pd.read_csv(input_csv)
    if "gene" not in table.columns:
        raise ValueError("Candidate table must contain a gene column.")
    genes = table["gene"].astype(str).str.upper()
    receptor = genes.isin(KNOWN_RECEPTORS) | genes.str.contains(r"^(?:GPR|ADGR|IL\d*R|CXCR|CCR|TNFRSF|TGFBR|FGFR|ERBB|PTPR|CD\d)", regex=True)
    transporter = genes.str.contains(r"^(?:SLC|ABCA|ABCB|ABCC|ATP\d|KCN|SCN|CACNA|CLCN)", regex=True)
    surface = genes.str.contains(r"^(?:CD\d|CLDN|ITGA|ITGB|EPCAM|MCAM|MUC\d)", regex=True)
    membrane = receptor | transporter | surface
    result = table.copy()
    result["receptor_like"] = receptor
    result["membrane_like"] = membrane
    result["transporter_like"] = transporter
    result["surface_like"] = surface
    result["annotation_reason"] = [
        ";".join(label for label, flag in (("receptor_like", r), ("transporter_like", t), ("surface_like", s)) if flag)
        for r, t, s in zip(receptor, transporter, surface)
    ]
    base_score = (
        pd.to_numeric(result["candidate_score"], errors="coerce").fillna(0.0)
        if "candidate_score" in result.columns else pd.Series(0.0, index=result.index)
    )
    result["priority_score"] = base_score + receptor.astype(float) + 0.5 * transporter.astype(float) + 0.5 * surface.astype(float)
    result = result[result["membrane_like"]].sort_values("priority_score", ascending=False)
    output = Path(output_csv); output.parent.mkdir(parents=True, exist_ok=True); result.to_csv(output, index=False)
    return result


def export_qubo_pool(input_csv: str | Path, output_csv: str | Path, max_genes: int = 40) -> list[str]:
    table = pd.read_csv(input_csv)
    if "gene" not in table.columns:
        raise ValueError("Input table must contain a gene column.")
    rank_column = "priority_score" if "priority_score" in table.columns else "candidate_score" if "candidate_score" in table.columns else None
    if rank_column:
        table = table.sort_values(rank_column, ascending=False)
    table = table.dropna(subset=["gene"]).copy()
    table["gene"] = table["gene"].astype(str)
    table = table.drop_duplicates(subset=["gene"], keep="first").head(max_genes)
    genes = table["gene"].tolist()
    if not genes:
        raise ValueError("No genes are available to export to the QUBO candidate pool.")
    output = Path(output_csv); output.parent.mkdir(parents=True, exist_ok=True)
    table.insert(1, "source_step", "A4" if "membrane_like" in table.columns else "A3")
    table["include_for_qubo"] = True
    preferred = [
        "gene", "source_step", "candidate_score", "direction", "receptor_like", "membrane_like",
        "transporter_like", "surface_like", "annotation_reason", "priority_score", "include_for_qubo",
    ]
    ordered = [column for column in preferred if column in table.columns]
    table[ordered + [column for column in table.columns if column not in ordered]].to_csv(output, index=False)
    return genes


def annotate_sequences(input_csv: str | Path, output_csv: str | Path) -> pd.DataFrame:
    table = pd.read_csv(input_csv)
    if "gene" not in table.columns:
        raise ValueError("Input table must contain a gene column.")
    result = table.copy()
    keys = result["gene"].astype(str).str.upper()
    result["normalized_gene"] = keys
    result["known_ligand"] = keys.isin(KNOWN_LIGANDS)
    result["known_receptor"] = keys.isin(KNOWN_RECEPTORS)
    for column in ("sequence_id", "sequence", "sequence_source", "sequence_review_status"):
        if column not in result:
            result[column] = ""
    output = Path(output_csv); output.parent.mkdir(parents=True, exist_ok=True); result.to_csv(output, index=False)
    return result


def build_ligrec_skeleton(input_csv: str | Path, output_csv: str | Path, max_pairs: int = 2500) -> pd.DataFrame:
    table = pd.read_csv(input_csv)
    if "gene" not in table.columns:
        raise ValueError("Input table must contain a gene column.")
    keys = table["gene"].astype(str).str.upper()
    ligands = list(dict.fromkeys(table.loc[keys.isin(KNOWN_LIGANDS), "gene"].astype(str)))
    receptors = list(dict.fromkeys(table.loc[keys.isin(KNOWN_RECEPTORS), "gene"].astype(str)))
    rows = []
    for ligand in ligands:
        for receptor in receptors:
            rows.append({"ligand_gene": ligand, "receptor_gene": receptor, "pair_status": "candidate_skeleton", "evidence_status": "not_reviewed"})
            if len(rows) >= max_pairs:
                break
        if len(rows) >= max_pairs:
            break
    result = pd.DataFrame(rows, columns=["ligand_gene", "receptor_gene", "pair_status", "evidence_status"])
    output = Path(output_csv); output.parent.mkdir(parents=True, exist_ok=True); result.to_csv(output, index=False)
    return result


def export_fasta_template(input_csv: str | Path, output_prefix: str | Path) -> tuple[Path, Path, int]:
    table = pd.read_csv(input_csv)
    prefix = Path(output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    template = prefix.with_suffix(".sequence_template.csv")
    fasta = prefix.with_suffix(".fasta")
    if "gene" not in table.columns:
        genes = pd.concat([table.get("ligand_gene", pd.Series(dtype=str)), table.get("receptor_gene", pd.Series(dtype=str))]).dropna().astype(str)
        table = pd.DataFrame({"gene": list(dict.fromkeys(genes))})
    result = table.copy()
    for column in ("sequence_id", "sequence", "sequence_source", "sequence_review_status"):
        if column not in result:
            result[column] = ""
    result.to_csv(template, index=False)
    records = []
    for _, row in result.iterrows():
        sequence = re.sub(r"[^A-Za-z*]", "", str(row.get("sequence", ""))).upper()
        if sequence and sequence != "NAN":
            identifier = str(row.get("sequence_id", "")).strip() or str(row["gene"])
            records.append(f">{identifier}|gene={row['gene']}\n{sequence}\n")
    fasta.write_text("".join(records), encoding="utf-8")
    return fasta, template, len(records)


def build_read_evidence_plan(input_csv: str | Path, output_csv: str | Path) -> pd.DataFrame:
    table = pd.read_csv(input_csv)
    if {"ligand_gene", "receptor_gene"}.issubset(table.columns):
        genes = list(dict.fromkeys(pd.concat([table["ligand_gene"], table["receptor_gene"]]).dropna().astype(str)))
    elif "gene" in table.columns:
        genes = list(dict.fromkeys(table["gene"].dropna().astype(str)))
    else:
        raise ValueError("Input must contain gene or ligand_gene/receptor_gene columns.")
    result = pd.DataFrame({
        "gene": genes, "expression_evidence": "review", "spatial_evidence": "review",
        "protein_localization_evidence": "review", "ligand_receptor_evidence": "review",
        "sequence_evidence": "review", "priority": "unassigned", "review_notes": "",
    })
    output = Path(output_csv); output.parent.mkdir(parents=True, exist_ok=True); result.to_csv(output, index=False)
    return result
