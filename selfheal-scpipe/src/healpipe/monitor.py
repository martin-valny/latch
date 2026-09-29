"""Post-step health checks.

Each check looks at a step's output and returns a CheckResult. Checks describe
*symptoms* ("no mitochondrial genes matched"), not causes. Working out the
cause is the agent's job.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
from anndata import AnnData

from .config import RunConfig
from .markers import MITO_PREFIX, marker_panel

BARCODE_RE = re.compile(r"^[ACGTN]{8,}(-\d+)?$")


@dataclass
class CheckResult:
    step: str
    name: str
    ok: bool
    message: str
    metrics: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {
            "step": self.step,
            "name": self.name,
            "ok": self.ok,
            "message": self.message,
            "metrics": self.metrics,
        }


def barcode_fraction(names) -> float:
    names = list(names[:500])
    if not names:
        return 0.0
    return sum(bool(BARCODE_RE.match(str(n))) for n in names) / len(names)


def _nonzero_values(X) -> np.ndarray:
    return X.data if sp.issparse(X) else X[X != 0]


def check_load(adata: AnnData, cfg: RunConfig) -> list[CheckResult]:
    vals = _nonzero_values(adata.X)
    n_nan = int(np.isnan(vals).sum())
    n_neg = int((vals < 0).sum())
    out = [
        CheckResult(
            "load",
            "finite_nonnegative",
            n_nan == 0 and n_neg == 0,
            "matrix is finite and non-negative"
            if n_nan == 0 and n_neg == 0
            else f"matrix has {n_nan} NaN and {n_neg} negative values",
            {"n_nan": n_nan, "n_negative": n_neg},
        )
    ]

    obs_bc = barcode_fraction(adata.obs_names)
    out.append(
        CheckResult(
            "load",
            "orientation",
            obs_bc >= 0.5,
            f"{obs_bc:.0%} of observation names look like cell barcodes",
            {"obs_barcode_fraction": round(obs_bc, 3), "shape": list(adata.shape)},
        )
    )

    finite = vals[np.isfinite(vals)]
    frac_int = float(np.mean(finite == np.round(finite))) if finite.size else 1.0
    vmax = float(finite.max()) if finite.size else 0.0
    if cfg.input_scale == "counts":
        ok = frac_int >= 0.999
        msg = (
            "values are integer counts"
            if ok
            else f"expected raw integer counts but {1 - frac_int:.1%} of non-zero values are non-integer"
        )
    else:
        ok = vmax < 20
        msg = f"log1p-scale input, max value {vmax:.2f}"
    out.append(
        CheckResult(
            "load",
            "value_scale",
            ok,
            msg,
            {"input_scale": cfg.input_scale, "frac_integer": round(frac_int, 4), "max": round(vmax, 3)},
        )
    )
    return out


def check_qc(adata: AnnData, cfg: RunConfig, n_cells_in: int) -> list[CheckResult]:
    prefix = MITO_PREFIX[cfg.species]
    n_mito = int(adata.uns["qc"]["n_mito_genes"])
    retained = adata.n_obs / max(n_cells_in, 1)
    return [
        CheckResult(
            "qc",
            "mito_genes_detected",
            n_mito > 0,
            f"{n_mito} mitochondrial genes matched prefix '{prefix}'",
            {"prefix": prefix, "n_mito_genes": n_mito},
        ),
        CheckResult(
            "qc",
            "cells_retained",
            retained >= 0.5,
            f"{adata.n_obs}/{n_cells_in} cells ({retained:.0%}) passed QC",
            {"n_in": n_cells_in, "n_out": int(adata.n_obs), "fraction": round(retained, 3)},
        ),
    ]


def check_normalize(adata: AnnData, cfg: RunConfig) -> list[CheckResult]:
    vals = _nonzero_values(adata.X)
    vmax = float(vals.max()) if vals.size else 0.0
    return [
        CheckResult(
            "normalize",
            "log_scale_range",
            0 < vmax < 20,
            f"max log-normalized value {vmax:.2f}",
            {"max": round(vmax, 3)},
        )
    ]


def check_annotate(adata: AnnData, cfg: RunConfig) -> list[CheckResult]:
    panel = marker_panel(cfg.species)
    wanted = {g for genes in panel.values() for g in genes}
    found = wanted & set(adata.var_names)
    overlap = len(found) / len(wanted)
    out = [
        CheckResult(
            "annotate",
            "marker_overlap",
            overlap >= 0.8,
            f"{len(found)}/{len(wanted)} panel markers present in the data",
            {"overlap": round(overlap, 3), "missing_examples": sorted(wanted - found)[:8]},
        )
    ]
    if "cell_type" in adata.obs:
        conf = float(adata.obs["annotation_confident"].mean())
        top = float(adata.obs["cell_type"].value_counts(normalize=True).iloc[0])
        out += [
            CheckResult(
                "annotate",
                "confident_fraction",
                conf >= 0.6,
                f"{conf:.0%} of cells confidently annotated",
                {"confident_fraction": round(conf, 3)},
            ),
            CheckResult(
                "annotate",
                "not_collapsed",
                top <= 0.9,
                f"largest cell type holds {top:.0%} of cells",
                {"largest_fraction": round(top, 3)},
            ),
        ]
    return out
