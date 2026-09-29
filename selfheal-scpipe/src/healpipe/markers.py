"""Canonical PBMC marker panel (textbook markers, e.g. the classic PBMC3k tutorial).

Keys match the reference labels shipped with the public PBMC3k dataset so the
demo can report agreement against them.
"""

from __future__ import annotations

HUMAN_PBMC_MARKERS: dict[str, list[str]] = {
    "CD4 T cells": ["CD3D", "CD3E", "IL7R", "LDHB"],
    "CD8 T cells": ["CD3D", "CD8A", "CD8B", "CCL5"],
    "NK cells": ["GNLY", "NKG7", "GZMB", "PRF1"],
    "B cells": ["MS4A1", "CD79A", "CD79B", "CD19"],
    "CD14+ Monocytes": ["CD14", "LYZ", "S100A8", "S100A9"],
    "FCGR3A+ Monocytes": ["FCGR3A", "MS4A7", "LST1", "IFITM3"],
    "Dendritic cells": ["FCER1A", "CLEC10A", "CST3"],
    "Megakaryocytes": ["PPBP", "PF4", "GNG11"],
}


def marker_panel(species: str) -> dict[str, list[str]]:
    if species == "human":
        return HUMAN_PBMC_MARKERS
    if species == "mouse":
        # Mouse symbols are title-case. This is an approximation of true
        # orthology, which is fine for a demo of the failure mode.
        return {
            ct: [g[0] + g[1:].lower() for g in genes]
            for ct, genes in HUMAN_PBMC_MARKERS.items()
        }
    raise ValueError(f"unknown species: {species}")


MITO_PREFIX = {"human": "MT-", "mouse": "mt-"}
