from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image
from scipy import sparse
from scipy.io import mmwrite


def write_flat_sample(
    root: Path,
    prefix: str,
    *,
    compressed: bool = False,
    include_images: bool = True,
    header_positions: bool = True,
    feature_count: int = 2,
    barcode_count: int = 3,
    matrix_shape: tuple[int, int] | None = None,
) -> dict[str, Path]:
    root.mkdir(parents=True, exist_ok=True)
    shape = matrix_shape or (feature_count, barcode_count)
    values = np.arange(1, shape[0] * shape[1] + 1, dtype=np.int32).reshape(shape)
    matrix = sparse.coo_matrix(values)
    matrix_path = root / f"{prefix}_matrix.mtx{'.gz' if compressed else ''}"
    if compressed:
        with gzip.open(matrix_path, "wb") as handle:
            mmwrite(handle, matrix)
    else:
        # A binary handle avoids platform-specific filename encoding behavior in SciPy.
        with matrix_path.open("wb") as handle:
            mmwrite(handle, matrix)

    barcodes = [f"BC{index + 1}" for index in range(barcode_count)]
    features = [f"GENE{index + 1}" for index in range(feature_count)]
    barcode_path = root / f"{prefix}_barcodes.tsv{'.gz' if compressed else ''}"
    feature_path = root / f"{prefix}_features.tsv{'.gz' if compressed else ''}"
    opener = gzip.open if compressed else open
    with opener(barcode_path, "wt", encoding="utf-8") as handle:
        handle.write("\n".join(barcodes) + "\n")
    with opener(feature_path, "wt", encoding="utf-8") as handle:
        handle.write("".join(f"ENSG{index + 1}\t{name}\tGene Expression\n" for index, name in enumerate(features)))

    positions = pd.DataFrame({
        "barcode": barcodes,
        "in_tissue": [1] * barcode_count,
        "array_row": list(range(barcode_count)),
        "array_col": list(reversed(range(barcode_count))),
        "pxl_row_in_fullres": [100 + index * 20 for index in range(barcode_count)],
        "pxl_col_in_fullres": [110 + index * 20 for index in range(barcode_count)],
    })
    position_path = root / f"{prefix}_{'tissue_positions' if header_positions else 'tissue_positions_list'}.csv{'.gz' if compressed else ''}"
    positions.to_csv(position_path, index=False, header=header_positions, compression="gzip" if compressed else None)
    scalefactor_path = root / f"{prefix}_scalefactors_json.json"
    scalefactor_path.write_text(json.dumps({
        "tissue_hires_scalef": 0.5,
        "tissue_lowres_scalef": 0.1,
        "spot_diameter_fullres": 80.0,
    }), encoding="utf-8")
    result = {
        "matrix": matrix_path,
        "barcodes": barcode_path,
        "features": feature_path,
        "positions": position_path,
        "scalefactors": scalefactor_path,
    }
    if include_images:
        hires = root / f"{prefix}_tissue_hires_image.png"
        lowres = root / f"{prefix}_tissue_lowres_image.png"
        Image.new("RGB", (4, 3), color=(20, 40, 60)).save(hires)
        Image.new("RGB", (2, 2), color=(70, 80, 90)).save(lowres)
        result.update({"hires": hires, "lowres": lowres})
    return result
