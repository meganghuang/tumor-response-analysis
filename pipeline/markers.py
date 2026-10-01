"""Data-driven marker genes per population."""
from __future__ import annotations

import pandas as pd
import scanpy as sc


def rank_markers(adata, groupby: str = "population", method: str = "wilcoxon") -> pd.DataFrame:
    """One-vs-rest ranked marker genes for every group, as a tidy table."""
    sub = adata[adata.obs[groupby].astype(str) != "Unresolved"].copy()
    sub.obs[groupby] = sub.obs[groupby].astype("category").cat.remove_unused_categories()
    sc.tl.rank_genes_groups(sub, groupby=groupby, method=method, pts=True,
                            key_added="markers", use_raw=False)
    df = sc.get.rank_genes_groups_df(sub, group=None, key="markers")
    df = df.rename(columns={
        "group": "population", "names": "gene", "logfoldchanges": "log2FC",
        "pvals": "pval", "pvals_adj": "pval_adj", "scores": "score",
        "pct_nz_group": "pct_in_group", "pct_nz_reference": "pct_in_rest",
    })
    df["rank"] = df.groupby("population").cumcount() + 1
    cols = ["population", "rank", "gene", "score", "log2FC", "pct_in_group",
            "pct_in_rest", "pval", "pval_adj"]
    return df[[c for c in cols if c in df.columns]]


def export_markers(markers: pd.DataFrame, alpha: float = 0.05, min_pct: float = 0.10
                   ) -> pd.DataFrame:
    """Significant, enriched, reasonably detected markers -- the publishable table.

    The full one-vs-rest table is every gene x every population (half a million
    rows, 64 MB) and is not a useful artifact. This keeps genes that are
    up-regulated in the population at BH q < ``alpha`` and detected in at least
    ``min_pct`` of its cells.
    """
    sel = markers[(markers["pval_adj"] < alpha) & (markers["log2FC"] > 0)
                  & (markers["pct_in_group"] >= min_pct)].copy()
    sel["rank"] = sel.groupby("population").cumcount() + 1
    return sel.reset_index(drop=True)


def top_markers(markers: pd.DataFrame, n: int = 5, min_pct: float = 0.25) -> dict[str, list[str]]:
    """Top-n significant markers per population, de-duplicated across populations."""
    sel = markers[(markers["pval_adj"] < 0.05) & (markers.get("pct_in_group", 1) >= min_pct)]
    out, seen = {}, set()
    for pop, grp in sel.groupby("population", observed=True):
        genes = [g for g in grp.sort_values("score", ascending=False)["gene"] if g not in seen]
        out[str(pop)] = genes[:n]
        seen.update(out[str(pop)])
    return out
