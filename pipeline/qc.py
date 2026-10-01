"""Per-cell quality control."""
from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc

MITO_PREFIX = "MT-"
RIBO_PREFIX = ("RPS", "RPL")


def annotate_qc(adata) -> None:
    """Attach n_genes, total TPM, and mitochondrial/ribosomal fractions to ``.obs``."""
    adata.var["mito"] = adata.var_names.str.startswith(MITO_PREFIX)
    adata.var["ribo"] = adata.var_names.str.startswith(RIBO_PREFIX)
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mito", "ribo"], percent_top=None, log1p=False, inplace=True
    )
    adata.obs["n_genes"] = adata.obs["n_genes_by_counts"]
    adata.obs["total_tpm"] = adata.obs["total_counts"]
    adata.obs["pct_mito"] = adata.obs["pct_counts_mito"]
    adata.obs["pct_ribo"] = adata.obs["pct_counts_ribo"]


def apply_qc(adata, cfg):
    """Filter cells, genes and under-sampled biopsies. Returns (adata, summary)."""
    q = cfg.qc
    annotate_qc(adata)
    n0, g0 = adata.shape

    keep = (
        (adata.obs["n_genes"] >= q["min_genes"])
        & (adata.obs["n_genes"] <= q["max_genes"])
        & (adata.obs["pct_mito"] <= q["max_pct_mito"])
    )
    steps = [
        ("input", n0),
        ("min_genes", int((adata.obs["n_genes"] >= q["min_genes"]).sum())),
        ("max_genes", int((adata.obs["n_genes"] <= q["max_genes"]).sum())),
        ("max_pct_mito", int((adata.obs["pct_mito"] <= q["max_pct_mito"]).sum())),
        ("all cell filters", int(keep.sum())),
    ]
    adata = adata[keep].copy()

    counts = adata.obs["sample"].value_counts()
    small = counts[counts < q["min_cells_per_sample"]].index
    adata = adata[~adata.obs["sample"].isin(small)].copy()
    steps.append((f"sample >= {q['min_cells_per_sample']} cells", adata.n_obs))

    sc.pp.filter_genes(adata, min_cells=q["min_cells_per_gene"])
    for col in ["sample", "patient", "timepoint", "response", "therapy"]:
        adata.obs[col] = adata.obs[col].cat.remove_unused_categories()

    summary = pd.DataFrame(steps, columns=["step", "n_cells"])
    summary["n_genes_retained"] = [g0] * (len(steps) - 1) + [adata.n_vars]
    return adata, summary
