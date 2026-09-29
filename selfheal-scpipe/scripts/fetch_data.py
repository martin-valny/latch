"""Download the public PBMC3k dataset and rebuild its raw integer count matrix.

Source: the example dataset shipped with CZI's cellxgene viewer (10x Genomics
PBMC3k, processed per the classic scanpy/Seurat tutorial). Its `.raw` holds
per-cell normalized log1p values. Dividing each cell's expm1 values by its
smallest non-zero value recovers the original UMI counts exactly, since a count
of 1 maps to that smallest value.

Usage:  python scripts/fetch_data.py [--out data/pbmc3k_counts.h5ad]
"""

from __future__ import annotations

import argparse
import tempfile
import urllib.request
from pathlib import Path

import anndata as ad
import numpy as np
import scipy.sparse as sp

URL = "https://raw.githubusercontent.com/chanzuckerberg/cellxgene/main/example-dataset/pbmc3k.h5ad"


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--out", default="data/pbmc3k_counts.h5ad")
    args = p.parse_args()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as tmp:
        src = Path(tmp) / "pbmc3k.h5ad"
        print(f"downloading {URL}")
        urllib.request.urlretrieve(URL, src)
        processed = ad.read_h5ad(src)

    raw = processed.raw.to_adata()
    X = raw.X.tocsr().astype(np.float64)
    X.data = np.expm1(X.data)
    per_cell_min = np.array([X.data[X.indptr[i] : X.indptr[i + 1]].min() for i in range(X.shape[0])])
    counts = sp.diags(1.0 / per_cell_min) @ X
    residual = np.abs(counts.data - np.round(counts.data)).max()
    assert residual < 1e-3, f"count reconstruction failed (max residual {residual})"
    counts.data = np.round(counts.data)
    counts = sp.csr_matrix(counts, dtype=np.float32)

    totals = np.asarray(counts.sum(axis=1)).ravel()
    assert np.allclose(totals, processed.obs["n_counts"].values), "totals do not match published n_counts"

    adata = ad.AnnData(X=counts, obs=processed.obs[[]].copy(), var=raw.var[[]].copy())
    adata.obs["reference_label"] = processed.obs["louvain"].astype(str).values
    adata.uns["source"] = {"url": URL, "description": "10x PBMC3k; counts rebuilt from normalized .raw"}
    adata.write_h5ad(out, compression="gzip")
    print(f"wrote {out}: {adata.n_obs} cells x {adata.n_vars} genes")


if __name__ == "__main__":
    main()
