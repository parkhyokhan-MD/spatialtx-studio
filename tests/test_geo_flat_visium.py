from __future__ import annotations

import gzip
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite

from spatialtx_desktop.importers.geo_flat_visium import (
    GeoFlatSample,
    convert_geo_flat_directory,
    convert_geo_flat_sample,
    scan_geo_flat_directory,
    validate_geo_flat_sample,
    write_comparative_manifest,
)
from spatialtx_desktop.importers.mex_to_h5ad import detect_mex_sample
from tests.geo_flat_fixtures import write_flat_sample


def hashes(root: Path) -> dict[str, str]:
    return {
        path.relative_to(root).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in root.rglob("*") if path.is_file()
    }


class GeoFlatVisiumTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_valid_sample_arbitrary_paths_conversion_provenance_and_source_safety(self) -> None:
        source = self.root / "임의 원본 폴더"
        output = self.root / "다른 출력 위치"
        prefix = "GSM9452684_sample_30_pre"
        write_flat_sample(source, prefix)
        before = hashes(source)
        samples = scan_geo_flat_directory(source.as_posix())
        self.assertEqual([sample.sample_prefix for sample in samples], [prefix])
        sample = samples[0]
        self.assertTrue(sample.valid, sample.errors)
        self.assertEqual(sample.geo_accession, "GSM9452684")
        self.assertEqual(sample.parsed_subject_id, "sample_30")
        self.assertEqual(sample.parsed_condition, "pre")
        converted, report = convert_geo_flat_sample(sample, output)
        self.assertTrue(report["valid"], report["errors"])
        self.assertEqual(before, hashes(source))
        self.assertEqual(set(before), {path.relative_to(source).as_posix() for path in source.iterdir()})
        adata = ad.read_h5ad(converted)
        self.assertEqual(adata.shape, (3, 2))
        self.assertIn("spatial", adata.obsm)
        provenance = adata.uns["spatialtx_import"]
        self.assertEqual(provenance["source_mode"], "geo_flat_visium")
        self.assertEqual(provenance["sample_prefix"], prefix)
        self.assertFalse(provenance["source_files_modified"])
        self.assertNotIn(str(source.parent), json.dumps(provenance, default=str))
        self.assertEqual(provenance["orientation_detected"], "genes_x_barcodes")
        self.assertEqual(provenance["orientation_action"], "transpose_to_barcodes_x_genes")
        image_refs = adata.uns["spatial"][prefix]["metadata"]["image_files_relative_to_h5ad_directory"]
        self.assertTrue(all(not Path(value).is_absolute() for value in image_refs))
        with self.assertRaises(FileExistsError):
            convert_geo_flat_sample(sample, output)

    def test_multiple_samples_full_prefix_grouping_metadata_and_deterministic_order(self) -> None:
        source = self.root / "flat"
        for prefix in (
            "GSM3_sample_3_post",
            "GSM2_sample_30_pre",
            "GSM1_sample_3_pre",
            "GSM4_sample_30_post",
        ):
            write_flat_sample(source, prefix, include_images=False)
        first = scan_geo_flat_directory(source)
        second = scan_geo_flat_directory(source)
        names = [sample.sample_prefix for sample in first]
        self.assertEqual(names, sorted(names, key=str.casefold))
        self.assertEqual(names, [sample.sample_prefix for sample in second])
        self.assertEqual(len(names), 4)
        self.assertIn("GSM1_sample_3_pre", names)
        self.assertIn("GSM3_sample_3_post", names)
        self.assertIn("GSM2_sample_30_pre", names)
        self.assertIn("GSM4_sample_30_post", names)
        self.assertEqual({sample.parsed_condition for sample in first}, {"pre", "post"})

    def test_compressed_and_uncompressed_components_and_header_layouts(self) -> None:
        source = self.root / "mixed components"
        write_flat_sample(source, "A_sample_1_pre", compressed=False, header_positions=True)
        write_flat_sample(source, "B_sample_1_post", compressed=True, header_positions=False)
        samples = scan_geo_flat_directory(source)
        self.assertEqual(len(samples), 2)
        self.assertTrue(all(sample.valid for sample in samples), [sample.errors for sample in samples])
        outputs = self.root / "outputs"
        result = convert_geo_flat_directory(source, outputs)
        self.assertEqual(set(result["converted_paths"]), {sample.sample_prefix for sample in samples})
        report_dir = result["report_dir"]
        self.assertTrue((report_dir / "geo_flat_inventory.csv").is_file())
        self.assertTrue((report_dir / "geo_flat_validation_report.csv").is_file())
        self.assertTrue((report_dir / "geo_flat_import_log.json").is_file())
        self.assertTrue((report_dir / "geo_flat_conversion_summary.csv").is_file())

    def test_repository_gzip_images_scalefactors_and_indexed_positions(self) -> None:
        source = self.root / "repository gzip"
        prefix = "GSM55_sample_55_pre"
        files = write_flat_sample(source, prefix, compressed=True, include_images=True)
        positions = pd.read_csv(files["positions"], compression="gzip")
        positions.to_csv(files["positions"], index=True, index_label="", compression="gzip")
        for key in ("scalefactors", "hires", "lowres"):
            path = files[key]
            target = path.with_name(path.name + ".gz")
            with path.open("rb") as source_handle, gzip.open(target, "wb") as target_handle:
                target_handle.write(source_handle.read())
            path.unlink()
        sample = scan_geo_flat_directory(source)[0]
        self.assertTrue(sample.valid, sample.errors)
        self.assertEqual(sample.coordinate_count, 3)
        self.assertTrue(sample.scalefactors_file.name.endswith(".json.gz"))
        self.assertEqual(len(sample.image_files), 2)
        output, report = convert_geo_flat_sample(sample, self.root / "repository output")
        self.assertTrue(report["valid"], report["errors"])
        self.assertTrue(output.is_file())

    def test_missing_required_components_are_explicit(self) -> None:
        for component in ("matrix", "barcodes", "features", "positions", "scalefactors"):
            with self.subTest(component=component):
                source = self.root / component
                files = write_flat_sample(source, f"X_sample_1_{component}")
                files[component].unlink()
                sample = scan_geo_flat_directory(source)[0]
                self.assertFalse(sample.valid)
                self.assertTrue(any("Missing required" in error for error in sample.errors), sample.errors)

    def test_missing_images_warns_but_converts_with_spatial_coordinates(self) -> None:
        source = self.root / "no images"
        write_flat_sample(source, "GSM8_sample_8_pre", include_images=False)
        sample = scan_geo_flat_directory(source)[0]
        self.assertTrue(sample.valid)
        self.assertEqual(sample.validation_status, "warning")
        self.assertTrue(any("image overlay is unavailable" in warning for warning in sample.warnings))
        output, _report = convert_geo_flat_sample(sample, self.root / "converted")
        adata = ad.read_h5ad(output)
        self.assertIn("spatial", adata.obsm)
        self.assertFalse(adata.uns["spatial"][sample.sample_prefix]["metadata"]["image_overlay_available"])

    def test_duplicate_matrix_and_position_candidates_are_rejected_without_guessing(self) -> None:
        source = self.root / "duplicates"
        prefix = "GSM9_sample_9_pre"
        files = write_flat_sample(source, prefix)
        with files["matrix"].open("rb") as source_handle, gzip.open(source / f"{prefix}_matrix.mtx.gz", "wb") as target:
            target.write(source_handle.read())
        positions = pd.read_csv(files["positions"])
        positions.to_csv(source / f"{prefix}_tissue_positions_list.csv", index=False, header=False)
        sample = scan_geo_flat_directory(source)[0]
        self.assertFalse(sample.valid)
        self.assertTrue(any("Duplicate matrix" in error for error in sample.errors), sample.errors)
        self.assertTrue(any("Duplicate positions" in error for error in sample.errors), sample.errors)

    def test_duplicate_barcode_feature_and_scalefactor_candidates_are_rejected(self) -> None:
        source = self.root / "more duplicates"
        prefix = "GSM10_sample_10_pre"
        files = write_flat_sample(source, prefix)
        with gzip.open(source / f"{prefix}_barcodes.tsv.gz", "wt", encoding="utf-8") as handle:
            handle.write(files["barcodes"].read_text(encoding="utf-8"))
        (source / f"{prefix}_genes.tsv").write_bytes(files["features"].read_bytes())
        with gzip.open(source / f"{prefix}_scalefactors_json.json.gz", "wt", encoding="utf-8") as handle:
            handle.write(files["scalefactors"].read_text(encoding="utf-8"))
        sample = scan_geo_flat_directory(source)[0]
        self.assertFalse(sample.valid)
        joined = "; ".join(sample.errors)
        self.assertIn("Duplicate barcodes", joined)
        self.assertIn("Duplicate features", joined)
        self.assertIn("Duplicate scalefactors", joined)

    def test_case_ambiguous_prefix_and_optional_auxiliary_inventory(self) -> None:
        source = self.root / "case prefixes"
        write_flat_sample(source, "Case_sample_1_pre", include_images=False)
        auxiliary = source / "Case_sample_1_pre_annotations.csv"
        auxiliary.write_text("barcode,label\nBC1,A\n", encoding="utf-8")
        sample = scan_geo_flat_directory(source)[0]
        self.assertIn(auxiliary.resolve(), sample.optional_files)
        ambiguous = GeoFlatSample("Case", source, prefix_variants={"Case", "case"})
        validate_geo_flat_sample(ambiguous)
        self.assertTrue(any("differing only by case" in error for error in ambiguous.errors), ambiguous.errors)

    def test_matrix_barcode_and_feature_mismatches_are_distinguished(self) -> None:
        barcode_source = self.root / "barcode mismatch"
        files = write_flat_sample(barcode_source, "A_sample_1_pre")
        files["barcodes"].write_text("BC1\nBC2\nBC3\nBC4\n", encoding="utf-8")
        sample = scan_geo_flat_directory(barcode_source)[0]
        self.assertFalse(sample.valid)
        self.assertTrue(any("barcode count" in error or "dimensions" in error for error in sample.errors), sample.errors)

        feature_source = self.root / "feature mismatch"
        files = write_flat_sample(feature_source, "B_sample_1_pre")
        files["features"].write_text("E1\tG1\nE2\tG2\nE3\tG3\n", encoding="utf-8")
        sample = scan_geo_flat_directory(feature_source)[0]
        self.assertFalse(sample.valid)
        self.assertTrue(any("feature count" in error or "dimensions" in error for error in sample.errors), sample.errors)

    def test_malformed_positions_scalefactors_gzip_and_text_are_rejected(self) -> None:
        source = self.root / "malformed"
        files = write_flat_sample(source, "A_sample_1_pre", compressed=True)
        files["positions"].write_bytes(b"not,a,valid,position\n")
        files["scalefactors"].write_text("{broken", encoding="utf-8")
        files["barcodes"].write_bytes(b"\x80\x81")
        sample = scan_geo_flat_directory(source)[0]
        self.assertFalse(sample.valid)
        joined = "; ".join(sample.errors)
        self.assertIn("Malformed tissue-position", joined)
        self.assertIn("Malformed scalefactor JSON", joined)
        self.assertIn("Malformed barcode", joined)

    def test_coordinate_mismatch_and_duplicate_barcodes_are_rejected(self) -> None:
        source = self.root / "coordinate mismatch"
        files = write_flat_sample(source, "A_sample_1_pre")
        positions = pd.read_csv(files["positions"])
        positions.loc[2, "barcode"] = "BC2"
        positions.to_csv(files["positions"], index=False)
        sample = scan_geo_flat_directory(source)[0]
        self.assertFalse(sample.valid)
        self.assertTrue(any("duplicate barcodes" in error for error in sample.errors), sample.errors)

        source2 = self.root / "missing coordinate"
        files2 = write_flat_sample(source2, "B_sample_1_pre")
        positions2 = pd.read_csv(files2["positions"]).iloc[:2]
        positions2.to_csv(files2["positions"], index=False)
        sample2 = scan_geo_flat_directory(source2)[0]
        self.assertFalse(sample2.valid)
        self.assertEqual(sample2.missing_coordinate_count, 1)

    def test_empty_inputs_and_ambiguous_square_orientation_fail(self) -> None:
        source = self.root / "empty"
        files = write_flat_sample(source, "A_sample_1_pre")
        files["barcodes"].write_text("", encoding="utf-8")
        sample = scan_geo_flat_directory(source)[0]
        self.assertFalse(sample.valid)
        self.assertTrue(any("Empty barcode" in error for error in sample.errors), sample.errors)

        source_features = self.root / "empty features"
        feature_files = write_flat_sample(source_features, "A2_sample_1_pre")
        feature_files["features"].write_text("", encoding="utf-8")
        feature_sample = scan_geo_flat_directory(source_features)[0]
        self.assertFalse(feature_sample.valid)
        self.assertTrue(any("Empty feature" in error for error in feature_sample.errors), feature_sample.errors)

        source2 = self.root / "ambiguous"
        write_flat_sample(source2, "B_sample_1_pre", feature_count=2, barcode_count=2, matrix_shape=(2, 2))
        sample2 = scan_geo_flat_directory(source2)[0]
        self.assertFalse(sample2.valid)
        self.assertEqual(sample2.orientation_detected, "ambiguous_equal_dimensions")

        source3 = self.root / "zero matrix"
        files3 = write_flat_sample(source3, "C_sample_1_pre")
        mmwrite(files3["matrix"], sparse.coo_matrix((2, 3), dtype=np.int32))
        sample3 = scan_geo_flat_directory(source3)[0]
        self.assertFalse(sample3.valid)
        self.assertTrue(any("empty" in error.lower() for error in sample3.errors))

    def test_recursive_scan_is_opt_in_and_source_location_does_not_change_data(self) -> None:
        source = self.root / "series"
        nested = source / "nested sample folder"
        prefix = "GSM20_sample_20_pre"
        write_flat_sample(nested, prefix)
        self.assertEqual(scan_geo_flat_directory(source, recursive=False), [])
        recursive = scan_geo_flat_directory(source, recursive=True)
        self.assertEqual([sample.sample_prefix for sample in recursive], [prefix])

        relocated = self.root / "다른 위치"
        relocated.mkdir()
        for path in nested.iterdir():
            (relocated / path.name).write_bytes(path.read_bytes())
        first, _ = convert_geo_flat_sample(recursive[0], self.root / "out1")
        second_sample = scan_geo_flat_directory(relocated)[0]
        second, _ = convert_geo_flat_sample(second_sample, self.root / "out2")
        a = ad.read_h5ad(first)
        b = ad.read_h5ad(second)
        np.testing.assert_array_equal(a.X.toarray(), b.X.toarray())
        np.testing.assert_allclose(a.obsm["spatial"], b.obsm["spatial"])

    def test_unconfirmed_manifest_does_not_auto_pair_and_confirmed_manifest_requires_approval(self) -> None:
        source = self.root / "manifest source"
        write_flat_sample(source, "GSM30_sample_30_pre", include_images=False)
        write_flat_sample(source, "GSM31_sample_30_post", include_images=False)
        result = convert_geo_flat_directory(source, self.root / "converted")
        unconfirmed_path = write_comparative_manifest(
            result["samples"], result["converted_paths"], self.root / "manifest_unconfirmed.csv"
        )
        unconfirmed = pd.read_csv(unconfirmed_path, keep_default_na=False)
        self.assertTrue(unconfirmed["pair_id"].eq("").all())
        self.assertTrue(unconfirmed["group"].eq("").all())
        self.assertTrue(unconfirmed["pairing_source"].eq("filename_inference_unconfirmed").all())

        confirmed_path = write_comparative_manifest(
            result["samples"], result["converted_paths"], self.root / "manifest_confirmed.csv", confirmed=True
        )
        confirmed = pd.read_csv(confirmed_path, keep_default_na=False)
        self.assertTrue(confirmed["pair_id"].eq("sample_30").all())
        self.assertEqual(set(confirmed["group"]), {"pre", "post"})
        self.assertTrue(confirmed["pairing_source"].eq("user_confirmed").all())

    def test_batch_output_collision_is_reported_and_standard_mex_importer_still_detects(self) -> None:
        source = self.root / "flat source"
        prefix = "GSM40_sample_40_pre"
        write_flat_sample(source, prefix, include_images=False)
        output = self.root / "output"
        first = convert_geo_flat_directory(source, output)
        self.assertEqual(first["summary"].iloc[0]["status"], "success")
        second = convert_geo_flat_directory(source, output)
        self.assertEqual(second["summary"].iloc[0]["status"], "failed_conversion")
        self.assertIn("already exists", second["summary"].iloc[0]["error_message"])

        canonical = self.root / "canonical_mex"
        canonical.mkdir()
        (canonical / "matrix.mtx").write_bytes((source / f"{prefix}_matrix.mtx").read_bytes())
        (canonical / "barcodes.tsv").write_bytes((source / f"{prefix}_barcodes.tsv").read_bytes())
        (canonical / "features.tsv").write_bytes((source / f"{prefix}_features.tsv").read_bytes())
        self.assertTrue(detect_mex_sample(canonical)["valid"])

    def test_output_inside_read_only_source_is_rejected_before_writing(self) -> None:
        source = self.root / "read only source boundary"
        write_flat_sample(source, "GSM60_sample_60_pre", include_images=False)
        before = hashes(source)
        with self.assertRaisesRegex(ValueError, "outside"):
            convert_geo_flat_directory(source, source / "converted")
        self.assertEqual(before, hashes(source))
        self.assertFalse((source / "converted").exists())

    def test_production_importer_has_no_machine_specific_path_literals(self) -> None:
        production_root = Path(__file__).resolve().parents[1] / "spatialtx_desktop"
        forbidden = ("GSE316402", "Desktop/", "Desktop\\", "C:\\Users\\", "/Users/", "/home/")
        violations = []
        for path in production_root.rglob("*.py"):
            text = path.read_text(encoding="utf-8")
            for token in forbidden:
                if token in text:
                    violations.append(f"{path.relative_to(production_root)}:{token}")
        self.assertEqual(violations, [])


if __name__ == "__main__":
    unittest.main()
