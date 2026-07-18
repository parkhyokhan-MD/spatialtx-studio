from __future__ import annotations

import logging
import traceback
from pathlib import Path

from .frame26 import run_frame26
from .io import ensure_output_dirs, load_config, load_h5ad, write_metrics, write_qc, write_selected_genes, write_yaml
from .metadata import DISCLAIMER, VERSION


def setup_logger(log_path: Path) -> logging.Logger:
    logger = logging.getLogger("spatialtx_studio")
    logger.handlers.clear()
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    file_handler = logging.FileHandler(log_path, encoding="utf-8")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    return logger


def run_one(input_path: str | Path, output_dir: str | Path, gene_mode: str = "fixed", config_path: str | Path | None = None, sample: str | None = None, analysis: str = "frame26") -> dict:
    paths = ensure_output_dirs(output_dir)
    log_path = paths["output"] / "run_log.txt"
    logger = setup_logger(log_path)
    try:
        config = load_config(None, config_path)
        config["input"]["path"] = str(input_path)
        config["output"]["path"] = str(output_dir)
        config["gene_mode"] = gene_mode
        config["analysis"] = analysis
        config["software_version"] = VERSION
        config["disclaimer"] = DISCLAIMER
        write_yaml(paths["output"] / "run_config.yaml", config)

        logger.info("Starting SpatialTX Studio run")
        logger.info("Input: %s", input_path)
        logger.info("Output: %s", output_dir)
        logger.info("Gene mode: %s", gene_mode)
        logger.info("Analysis: %s", analysis)

        if analysis == "frame26":
            sample_name = sample or Path(input_path).stem
            summary, _adata = run_frame26(input_path, sample_name, paths["output"], mode=gene_mode, config=config)
            # run_frame26 writes its own detailed run_log; append top-level completion line.
            with log_path.open("a", encoding="utf-8") as f:
                f.write("analysis=frame26\n")
                f.write("status=ok\n")
            return summary.iloc[0].to_dict()
        if analysis != "istz":
            raise ValueError(f"Unsupported analysis: {analysis}")

        # Legacy ISTZ dependencies are loaded only when that explicitly
        # requested analysis is used. The canonical FRAME2.6 CLI and --help do
        # not require Scanpy at import time.
        from .distance_calibration import resolve_distance_config
        from .gene_program import select_gene_programs
        from .interface_detection import detect_interface
        from .plotting import save_interface_map, save_transition_zone_map
        from .preprocess import preprocess_adata
        from .qc import build_qc_report
        from .spatial_fields import build_fields
        from .transition_metrics import compute_metrics
        from .transition_zone import compute_transition_zone

        adata = load_h5ad(input_path)
        adata = preprocess_adata(
            adata,
            min_genes_per_spot=int(config["preprocessing"]["min_genes_per_spot"]),
            gene_min_spots=int(config["preprocessing"]["gene_min_spots"]),
            use_in_tissue=bool(config["preprocessing"]["use_in_tissue"]),
            target_sum=float(config["preprocessing"]["target_sum"]),
        )
        selected_c, selected_b, selection_note = select_gene_programs(adata, config, gene_mode)
        fields = build_fields(adata, selected_c, selected_b, config)
        interface = detect_interface(fields, config)
        distance_config = resolve_distance_config(config)
        transition = compute_transition_zone(fields, interface, config)
        sample_name = sample or Path(input_path).stem
        metrics = compute_metrics(sample_name, gene_mode, adata, selected_c, selected_b, fields, interface, transition)

        qc_rows = build_qc_report(metrics, selection_note, fields["coord_source"], distance_config["note"])
        write_metrics(paths["output"] / "metrics.csv", metrics)
        write_qc(paths["output"] / "qc_report.csv", qc_rows)
        write_selected_genes(paths["output"] / "selected_genes.csv", selected_c, selected_b)
        save_interface_map(paths["figures"] / "interface_map.png", fields, interface)
        save_transition_zone_map(paths["figures"] / "transition_zone_map.png", fields, interface, transition)
        logger.info("Run completed successfully")
        return metrics
    except Exception as exc:
        logger.error("Run failed: %s", exc)
        logger.error(traceback.format_exc())
        raise RuntimeError(f"SpatialTX Studio run failed. See log: {log_path}") from exc
