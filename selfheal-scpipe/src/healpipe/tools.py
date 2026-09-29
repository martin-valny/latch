"""The agent's tool surface. Every tool is deterministic Python.

The LLM chooses which tool to call and with what arguments, and the tools do
the actual work. Guardrails are enforced here, in code, not in the prompt:

* only allow-listed config keys and values can be changed (QC thresholds cannot)
* relaunches are capped
* a relaunch must follow a config change and must start early enough to
  recompute every step that change affects, so no blind retries
* inspection tools are read-only
"""

from __future__ import annotations

from dataclasses import dataclass, field

import anndata as ad
import numpy as np
import scipy.sparse as sp

from .config import AGENT_EDITABLE, STEPS
from .markers import MITO_PREFIX
from .monitor import barcode_fraction
from .pipeline import Pipeline
from .trace import Trace


class ToolError(Exception):
    pass


@dataclass
class Toolbox:
    pipeline: Pipeline
    trace: Trace
    max_relaunches: int = 3
    relaunches: int = 0
    dirty_from: str | None = None  # earliest step invalidated by config edits
    changes: list[dict] = field(default_factory=list)
    escalation: dict | None = None

    # ------------------------------------------------------------ dispatch

    def call(self, name: str, args: dict) -> dict:
        self.trace.add("tool_call", tool=name, input=args)
        fn = getattr(self, f"tool_{name}", None)
        try:
            if fn is None:
                raise ToolError(f"unknown tool '{name}'")
            params = {k: v for k, v in args.items() if k != "rationale"}
            result = fn(**params)
            ok = True
        except (ToolError, TypeError, ValueError) as e:
            result, ok = {"error": str(e)}, False
        self.trace.add("tool_result", tool=name, ok=ok, result=result)
        return result

    @property
    def done(self) -> bool:
        passed = self.pipeline.last_result is not None and self.pipeline.last_result.status == "passed"
        return self.escalation is not None or passed

    # ------------------------------------------------------------ read-only

    def tool_get_run_status(self) -> dict:
        r = self.pipeline.last_result
        return {
            "config": self.pipeline.cfg.to_dict(),
            "agent_editable": {k: sorted(v) for k, (v, _) in AGENT_EDITABLE.items()},
            "last_run": r.to_dict() if r else None,
            "relaunches_used": self.relaunches,
            "relaunches_left": self.max_relaunches - self.relaunches,
        }

    def tool_inspect_matrix(self) -> dict:
        """Facts about the raw input file exactly as stored on disk, before any orientation fix."""
        a = ad.read_h5ad(self.pipeline.input_path)
        X = a.X if sp.issparse(a.X) else sp.csr_matrix(a.X)
        X = X.tocsr()
        v = X.data
        finite = v[np.isfinite(v)]
        out = {
            "shape_on_disk": list(a.shape),
            "obs_name_examples": list(map(str, a.obs_names[:4])),
            "var_name_examples": list(map(str, a.var_names[:4])),
            "obs_barcode_like_fraction": round(barcode_fraction(a.obs_names), 3),
            "var_barcode_like_fraction": round(barcode_fraction(a.var_names), 3),
            "density": round(X.nnz / (X.shape[0] * X.shape[1]), 4),
            "n_nan": int(np.isnan(v).sum()),
            "n_negative": int((v < 0).sum()),
            "min_nonzero": float(finite.min()) if finite.size else None,
            "max": float(finite.max()) if finite.size else None,
            "frac_nonzero_integer": round(float(np.mean(finite == np.round(finite))), 4) if finite.size else None,
        }
        # If values look log-scaled: after expm1, are per-row totals ~constant?
        # A constant total means the rows were library-size normalized before the log.
        cell_axis = 1 if out["var_barcode_like_fraction"] > out["obs_barcode_like_fraction"] else 0
        if finite.size and out["frac_nonzero_integer"] < 0.999 and out["n_negative"] == 0:
            E = X.copy()
            E.data = np.expm1(E.data)
            sums = np.asarray(E.sum(axis=1 - cell_axis)).ravel()
            sums = sums[sums > 0]
            out["expm1_per_cell_total_median"] = round(float(np.median(sums)), 1)
            out["expm1_per_cell_total_cv"] = round(float(sums.std() / sums.mean()), 4)
        return out

    def tool_inspect_gene_names(self) -> dict:
        """Gene nomenclature on the gene axis, taking the current orientation config into account."""
        a = ad.read_h5ad(self.pipeline.input_path, backed="r")
        genes = a.obs_names if self.pipeline.cfg.orientation == "genes_x_cells" else a.var_names
        genes = [str(g) for g in genes]
        letters = [g for g in genes if any(c.isalpha() for c in g)]
        upper = sum(g == g.upper() for g in letters) / max(len(letters), 1)
        title = sum(g[:1].isupper() and g[1:] == g[1:].lower() for g in letters) / max(len(letters), 1)
        n_human_mt = sum(g.startswith(MITO_PREFIX["human"]) for g in genes)
        n_mouse_mt = sum(g.startswith(MITO_PREFIX["mouse"]) for g in genes)
        return {
            "n_genes": len(genes),
            "examples": genes[:3] + [g for g in genes if g.upper().startswith("MT-")][:3],
            "fraction_all_uppercase": round(upper, 3),
            "fraction_titlecase": round(title, 3),
            "n_prefixed_MT-": n_human_mt,
            "n_prefixed_mt-": n_mouse_mt,
            "note": "Human symbols are conventionally uppercase (MT-CO1). Mouse symbols are title-case (mt-Co1).",
        }

    def tool_inspect_qc_distribution(self) -> dict:
        """Per-cell depth and complexity of the loaded matrix, versus the (fixed) QC thresholds."""
        a = self.pipeline.cache.get("load")
        if a is None:
            raise ToolError("the load step has not passed yet, so there is no loaded matrix to inspect")
        X = a.X
        if self.pipeline.cfg.input_scale == "log1p":
            X = X.copy()
            X.data = np.expm1(X.data)
        n_genes = np.asarray((a.X > 0).sum(axis=1)).ravel()
        total = np.asarray(X.sum(axis=1)).ravel()
        q = lambda x: {p: round(float(np.percentile(x, p)), 1) for p in (5, 25, 50, 75, 95)}
        cfg = self.pipeline.cfg
        return {
            "n_cells": int(a.n_obs),
            "genes_per_cell_percentiles": q(n_genes),
            "umis_per_cell_percentiles": q(total),
            "min_genes_threshold": cfg.min_genes,
            "fraction_cells_passing_min_genes": round(float((n_genes >= cfg.min_genes).mean()), 3),
            "note": "QC thresholds are fixed by the analysis owner and are not agent-editable.",
        }

    # ------------------------------------------------------------ actions

    def tool_set_config(self, key: str, value: str) -> dict:
        if key not in AGENT_EDITABLE:
            raise ToolError(
                f"'{key}' is not agent-editable. Editable keys: {sorted(AGENT_EDITABLE)}. "
                "Thresholds are owned by humans. If a threshold looks like the problem, escalate."
            )
        allowed, first_step = AGENT_EDITABLE[key]
        if value not in allowed:
            raise ToolError(f"invalid value for {key}: {value!r}; allowed: {sorted(allowed)}")
        old = getattr(self.pipeline.cfg, key)
        if old == value:
            raise ToolError(f"{key} is already {value!r}")
        self.pipeline.cfg = self.pipeline.cfg.with_(**{key: value})
        self.changes.append({"key": key, "old": old, "new": value})
        if self.dirty_from is None or STEPS.index(first_step) < STEPS.index(self.dirty_from):
            self.dirty_from = first_step
        return {"changed": {key: [old, value]}, "must_relaunch_from_or_before": self.dirty_from}

    def tool_relaunch(self, from_step: str) -> dict:
        if from_step not in STEPS:
            raise ToolError(f"unknown step {from_step!r}; steps: {STEPS}")
        if self.relaunches >= self.max_relaunches:
            raise ToolError(f"relaunch budget exhausted ({self.max_relaunches}). Escalate.")
        if self.dirty_from is None:
            raise ToolError("no config change since the last run. Relaunching would reproduce the same failure.")
        if STEPS.index(from_step) > STEPS.index(self.dirty_from):
            raise ToolError(f"config changes affect '{self.dirty_from}'. Relaunch from '{self.dirty_from}' or earlier.")
        self.relaunches += 1
        result = self.pipeline.run(from_step)
        self.dirty_from = None
        d = result.to_dict()
        self.trace.add("pipeline_run", from_step=from_step, **d)
        d["relaunches_left"] = self.max_relaunches - self.relaunches
        return d

    def tool_escalate(self, reason: str, evidence: str = "") -> dict:
        self.escalation = {"reason": reason, "evidence": evidence}
        return {"escalated": True, "message": "A human has been notified. Stop here."}


def _prop(desc: str, enum: list[str] | None = None) -> dict:
    p = {"type": "string", "description": desc}
    if enum:
        p["enum"] = enum
    return p


_RATIONALE = _prop("One sentence: the hypothesis this call tests, or why you are taking this action.")


def _schema(props: dict | None = None) -> dict:
    props = {"rationale": _RATIONALE, **(props or {})}
    return {"type": "object", "properties": props, "required": list(props), "additionalProperties": False}


TOOL_SPECS: list[dict] = [
    {
        "name": "get_run_status",
        "description": "Current run config, which keys you may edit, the last run's failed and passed checks, and the remaining relaunch budget.",
        "input_schema": _schema(),
    },
    {
        "name": "inspect_matrix",
        "description": "Read-only facts about the input matrix as stored on disk: shape, axis name examples, barcode-likeness of each axis, value range, integer fraction, negatives/NaNs, and for log-looking data whether per-cell expm1 totals are constant (which suggests it was normalized before the log).",
        "input_schema": _schema(),
    },
    {
        "name": "inspect_gene_names",
        "description": "Read-only gene nomenclature summary: case conventions and counts of human-style (MT-) versus mouse-style (mt-) mitochondrial genes.",
        "input_schema": _schema(),
    },
    {
        "name": "inspect_qc_distribution",
        "description": "Read-only per-cell depth and complexity distributions against the fixed QC thresholds. Needs the load step to have passed.",
        "input_schema": _schema(),
    },
    {
        "name": "set_config",
        "description": "Change one agent-editable run-config key (species, input_scale, orientation). Other keys, including all QC thresholds, are rejected.",
        "input_schema": _schema(
            {
                "key": _prop("Config key.", sorted(AGENT_EDITABLE)),
                "value": _prop("New value.", sorted({v for vs, _ in AGENT_EDITABLE.values() for v in vs})),
            }
        ),
    },
    {
        "name": "relaunch",
        "description": "Re-run the pipeline from a step, after a config change. Must start at or before the earliest step the change affects. Budget-limited.",
        "input_schema": _schema({"from_step": _prop("Step to resume from.", STEPS)}),
    },
    {
        "name": "escalate",
        "description": "Stop and hand off to a human. Use it when the evidence points to a problem with the data itself, when the only fix would need a non-editable setting, or when you are not confident.",
        "input_schema": _schema(
            {
                "reason": _prop("Short diagnosis for the human."),
                "evidence": _prop("Key numbers that support the diagnosis."),
            }
        ),
    },
]

