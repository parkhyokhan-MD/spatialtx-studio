from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


COMPARATIVE_GENES = ["C1", "C2", "S1", "S2", "CA9", "VEGFA", "PECAM1", "VWF"]
C_GENES = ["C1", "C2"]
S_GENES = ["S1", "S2"]


def write_comparative_h5ad(
    path: Path,
    *,
    pattern: str = "localized",
    include_context: bool = True,
    spatial: bool = True,
    seed: int = 1,
) -> Path:
    rng = np.random.default_rng(seed)
    side = 6
    coords = np.asarray([[x, y] for y in range(side) for x in range(side)], dtype=float)
    x, y = coords[:, 0], coords[:, 1]
    matrix = rng.poisson(2.0, size=(len(coords), len(COMPARATIVE_GENES))).astype(float)
    if pattern == "localized":
        matrix[:, 0] += np.where(x <= 3, 9, 1)
        matrix[:, 1] += np.where(x <= 2, 7, 1)
        matrix[:, 2] += np.where(x >= 2, 9, 1)
        matrix[:, 3] += np.where(x >= 3, 7, 1)
    elif pattern == "diffuse":
        checker = ((x.astype(int) + y.astype(int)) % 2) == 0
        matrix[:, 0] += np.where(checker, 10, 1)
        matrix[:, 1] += np.where((x.astype(int) % 3) == 0, 8, 1)
        matrix[:, 2] += np.where(~checker, 10, 1)
        matrix[:, 3] += np.where((y.astype(int) % 3) == 0, 8, 1)
    elif pattern == "flat":
        matrix[:, :4] += 4
    else:
        raise ValueError(pattern)
    if include_context:
        matrix[:, 4] += np.where(y >= 3, 6, 0)
        matrix[:, 5] += np.where(y >= 2, 4, 0)
        matrix[:, 6] += np.where(x >= 3, 5, 0)
        matrix[:, 7] += np.where(x >= 4, 7, 0)
    else:
        matrix = matrix[:, :4]
        genes = COMPARATIVE_GENES[:4]
    genes = COMPARATIVE_GENES if include_context else COMPARATIVE_GENES[:4]
    adata = ad.AnnData(matrix, var=pd.DataFrame(index=genes))
    if spatial:
        adata.obsm["spatial"] = coords
    path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(path)
    return path
