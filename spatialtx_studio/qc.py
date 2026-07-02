from __future__ import annotations


def build_qc_report(metrics: dict, selection_note: str, coord_source: str, distance_note: str) -> list[dict]:
    rows = [
        {"check": "run_status", "status": metrics.get("qc_status", "unknown"), "detail": ""},
        {"check": "gene_selection", "status": "ok", "detail": selection_note},
        {"check": "spatial_coordinates", "status": "ok", "detail": coord_source},
        {"check": "interface_detection", "status": "ok" if metrics["interface_detected"] else "warning", "detail": f"interface_fraction={metrics['interface_fraction']}"},
        {"check": "distance_mode", "status": "ok", "detail": distance_note},
        {"check": "research_use_only", "status": "notice", "detail": "Not intended for clinical diagnosis or treatment decision-making."},
    ]
    return rows
