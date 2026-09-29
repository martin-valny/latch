"""A small synthetic PBMC-like count matrix, so tests run offline in seconds."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from healpipe.markers import HUMAN_PBMC_MARKERS


@pytest.fixture(scope="session")
def clean_h5ad(tmp_path_factory):
    rng = np.random.default_rng(0)
    types = list(HUMAN_PBMC_MARKERS)
    markers = sorted({g for gs in HUMAN_PBMC_MARKERS.values() for g in gs})
    mito = [f"MT-GENE{i}" for i in range(13)]
    filler = [f"GENE{i}" for i in range(700)]
    genes = markers + mito + filler

    n_cells = 480
    labels = np.array(types)[np.arange(n_cells) % len(types)]
    lam = np.full((n_cells, len(genes)), 0.6)
    lam[:, len(markers) : len(markers) + len(mito)] = 0.4
    lam[:, : len(markers)] = 0.1
    col = {g: i for i, g in enumerate(genes)}
    for i, ct in enumerate(labels):
        for g in HUMAN_PBMC_MARKERS[ct]:
            lam[i, col[g]] = 6.0
    X = sp.csr_matrix(rng.poisson(lam).astype(np.float32))

    barcodes = ["".join(rng.choice(list("ACGT"), 14)) + "-1" for _ in range(n_cells)]
    a = ad.AnnData(X=X, obs=pd.DataFrame({"reference_label": labels}, index=barcodes), var=pd.DataFrame(index=genes))
    path = tmp_path_factory.mktemp("data") / "synthetic_counts.h5ad"
    a.write_h5ad(path)
    return path
