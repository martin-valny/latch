"""Run configuration and the subset of it the agent is allowed to change."""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace

STEPS = ["load", "qc", "normalize", "annotate"]


@dataclass(frozen=True)
class RunConfig:
    # Sample-sheet style metadata. These are the things that are plausibly
    # *wrong* in a real run, so the agent may correct them.
    species: str = "human"  # "human" | "mouse"
    input_scale: str = "counts"  # "counts" | "log1p"
    orientation: str = "cells_x_genes"  # "cells_x_genes" | "genes_x_cells"

    # Analysis thresholds. These encode scientific judgement, so the agent
    # may NOT change them. Loosening QC to make a failing run "pass" is the
    # classic bad automated fix.
    min_genes: int = 200
    max_pct_mito: float = 5.0
    min_cells_per_gene: int = 3
    target_sum: float = 1e4

    def with_(self, **kw) -> "RunConfig":
        return replace(self, **kw)

    def to_dict(self) -> dict:
        return asdict(self)


# key -> (allowed values, first pipeline step whose output depends on it)
AGENT_EDITABLE: dict[str, tuple[set[str], str]] = {
    "species": ({"human", "mouse"}, "qc"),
    "input_scale": ({"counts", "log1p"}, "load"),
    "orientation": ({"cells_x_genes", "genes_x_cells"}, "load"),
}
