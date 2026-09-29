"""A small, linear single-cell pipeline: load -> qc -> normalize -> annotate.

Each step's output is cached, so a relaunch can resume from any step. After
each step the monitor runs, and the first failing check stops the run.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from anndata import AnnData

from . import monitor
from .config import STEPS, RunConfig
from .markers import MITO_PREFIX, marker_panel

warnings.filterwarnings("ignore", category=FutureWarning)


@dataclass
class RunResult:
    status: str  # "passed" | "failed"
    failed_step: str | None
    checks: list[monitor.CheckResult] = field(default_factory=list)
    steps_run: list[str] = field(default_factory=list)

    @property
    def failed_checks(self) -> list[monitor.CheckResult]:
        return [c for c in self.checks if not c.ok]

    def to_dict(self) -> dict:
        return {
            "status": self.status,
            "failed_step": self.failed_step,
            "steps_run": self.steps_run,
            "failed_checks": [c.to_dict() for c in self.failed_checks],
            "passed_checks": [f"{c.step}.{c.name}" for c in self.checks if c.ok],
        }


# ---------------------------------------------------------------- steps


def step_load(path: Path, cfg: RunConfig) -> AnnData:
    adata = ad.read_h5ad(path)
    if cfg.orientation == "genes_x_cells":
        adata = adata.T.copy()
    X = adata.X
    adata.X = (X if sp.issparse(X) else sp.csr_matrix(X)).astype(np.float32).tocsr()
    return adata


def step_qc(adata: AnnData, cfg: RunConfig) -> AnnData:
    adata = adata.copy()
    prefix = MITO_PREFIX[cfg.species]
    adata.var["mito"] = adata.var_names.str.startswith(prefix)
    counts = adata.X
    if cfg.input_scale == "log1p":
        # Mito % is a ratio, so it is unaffected by normalization only if we
        # undo the log. Library-size normalization cancels out in the ratio.
        counts = counts.copy()
        counts.data = np.expm1(counts.data)
    total = np.asarray(counts.sum(axis=1)).ravel()
    mito = np.asarray(counts[:, adata.var["mito"].values].sum(axis=1)).ravel()
    adata.obs["n_genes"] = np.asarray((adata.X > 0).sum(axis=1)).ravel()
    adata.obs["pct_mito"] = np.divide(100 * mito, total, out=np.zeros_like(total), where=total > 0)

    keep = (adata.obs["n_genes"] >= cfg.min_genes) & (adata.obs["pct_mito"] < cfg.max_pct_mito)
    adata = adata[keep.values].copy()
    if adata.n_obs:
        sc.pp.filter_genes(adata, min_cells=cfg.min_cells_per_gene)
    adata.uns["qc"] = {"n_mito_genes": int(adata.var.get("mito", pd.Series(dtype=bool)).sum())}
    return adata


def step_normalize(adata: AnnData, cfg: RunConfig) -> AnnData:
    adata = adata.copy()
    if cfg.input_scale == "counts":
        adata.layers["counts"] = adata.X.copy()
        sc.pp.normalize_total(adata, target_sum=cfg.target_sum)
        sc.pp.log1p(adata)
    return adata


def step_annotate(adata: AnnData, cfg: RunConfig) -> AnnData:
    """Marker-score annotation: z-score each marker, average per cell type, take argmax."""
    adata = adata.copy()
    panel = {
        ct: [g for g in genes if g in adata.var_names]
        for ct, genes in marker_panel(cfg.species).items()
    }
    panel = {ct: g for ct, g in panel.items() if g}
    if not panel or adata.n_obs == 0:
        return adata

    genes = sorted({g for gs in panel.values() for g in gs})
    X = adata[:, genes].X
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    sd = X.std(axis=0)
    Z = np.clip((X - X.mean(axis=0)) / np.where(sd > 0, sd, 1), -10, 10)
    col = {g: i for i, g in enumerate(genes)}
    types = list(panel)
    scores = np.column_stack([Z[:, [col[g] for g in panel[ct]]].mean(axis=1) for ct in types])

    order = np.argsort(scores, axis=1)
    best = scores[np.arange(len(scores)), order[:, -1]]
    second = scores[np.arange(len(scores)), order[:, -2]] if len(types) > 1 else best * 0
    adata.obs["cell_type"] = pd.Categorical([types[i] for i in order[:, -1]])
    adata.obs["annotation_margin"] = best - second
    adata.obs["annotation_confident"] = (best > 0.5) & (best - second > 0.25)
    return adata


# ---------------------------------------------------------------- runner


class Pipeline:
    def __init__(self, input_path: Path | str, cfg: RunConfig):
        self.input_path = Path(input_path)
        self.cfg = cfg
        self.cache: dict[str, AnnData] = {}
        self.last_result: RunResult | None = None
        self.n_runs = 0

    def run(self, from_step: str = "load") -> RunResult:
        start = STEPS.index(from_step)
        if start > 0 and STEPS[start - 1] not in self.cache:
            raise ValueError(f"cannot resume from '{from_step}': no cached output for '{STEPS[start - 1]}'")
        # Anything downstream of the resume point is stale.
        for s in STEPS[start:]:
            self.cache.pop(s, None)

        self.n_runs += 1
        result = RunResult(status="passed", failed_step=None)
        for step in STEPS[start:]:
            prev = self.cache.get(STEPS[STEPS.index(step) - 1]) if step != "load" else None
            if step == "load":
                out = step_load(self.input_path, self.cfg)
                checks = monitor.check_load(out, self.cfg)
            elif step == "qc":
                out = step_qc(prev, self.cfg)
                checks = monitor.check_qc(out, self.cfg, prev.n_obs)
            elif step == "normalize":
                out = step_normalize(prev, self.cfg)
                checks = monitor.check_normalize(out, self.cfg)
            else:
                out = step_annotate(prev, self.cfg)
                checks = monitor.check_annotate(out, self.cfg)

            result.steps_run.append(step)
            result.checks.extend(checks)
            if not all(c.ok for c in checks):
                result.status, result.failed_step = "failed", step
                break
            self.cache[step] = out

        self.last_result = result
        return result

    def reference_agreement(self) -> float | None:
        """Agreement with the dataset's published labels. For reporting only; the monitor never sees it."""
        a = self.cache.get("annotate")
        if a is None or "reference_label" not in a.obs or "cell_type" not in a.obs:
            return None
        return float((a.obs["reference_label"].astype(str) == a.obs["cell_type"].astype(str)).mean())
