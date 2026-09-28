import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

# Markers for the CD45+ immune populations described in Sade-Feldman et al. 2018.
# Genes missing from the dataset are skipped (e.g. JCHAIN is named IGJ in this annotation).
# Interferon/activation-induced genes (IRF7, LAMP3) are avoided: activated T cells express them.
MARKER_GENES = {
    "B cells": ["CD19", "MS4A1", "CD79A", "CD79B"],
    "Plasma cells": ["IGJ", "JCHAIN", "MZB1", "SDC1"],
    "Monocytes-Macrophages": ["CD14", "CD68", "LYZ", "C1QA", "C1QB", "CSF1R"],
    "Dendritic cells": ["FCER1A", "CLEC10A", "CD1C"],
    "Plasmacytoid DCs": ["LILRA4", "CLEC4C", "IL3RA"],
    # Effector vs exhausted CD8 states are both granzyme-high, so granzymes and CD8A (shared)
    # are left out and each state is defined by what sets it apart.
    "CD8 effector": ["GZMH", "GNLY", "FGFBP2", "CX3CR1", "KLRG1"],
    "CD8 exhausted": ["PDCD1", "HAVCR2", "LAG3", "ENTPD1", "TOX", "CXCL13"],
    "Memory-naive T": ["TCF7", "IL7R", "CCR7", "SELL", "LEF1"],
    "Tregs": ["FOXP3", "IL2RA", "CTLA4"],
    # KLRD1 is left out: CD8 T cells express it too.
    "NK cells": ["NCAM1", "KLRF1", "SH2D1B"],
    "Cycling": ["MKI67", "TOP2A", "TYMS"],
}

# Populations that cannot contain T cells; used to find CD8 T cells in any other cluster.
NON_T_CELL_TYPES = {"B cells", "Plasma cells", "Monocytes-Macrophages", "Dendritic cells",
                    "Plasmacytoid DCs", "NK cells"}
UNRESOLVED = "Unresolved"


def cluster(adata: ad.AnnData, resolution: float) -> ad.AnnData:
    sc.tl.leiden(adata, resolution=resolution, key_added="leiden", flavor="igraph", n_iterations=2)
    return adata


def available_markers(adata: ad.AnnData, markers: dict[str, list[str]] = MARKER_GENES) -> dict[str, list[str]]:
    present = {k: [g for g in v if g in adata.raw.var_names] for k, v in markers.items()}
    return {k: v for k, v in present.items() if v}


def score_cell_types(
    adata: ad.AnnData,
    markers: dict[str, list[str]] = MARKER_GENES,
    min_score: float = 0.2,
    min_margin: float = 0.1,
) -> ad.AnnData:
    """Score each cell for each marker set and label clusters by their top-scoring cell type.

    A cluster whose best mean score is below `min_score`, or beats the runner-up by less than
    `min_margin`, is labelled "Unresolved" instead of being forced into a weak match.
    """
    score_cols = []
    for cell_type, genes in available_markers(adata, markers).items():
        col = f"score_{cell_type}"
        sc.tl.score_genes(adata, genes, score_name=col, use_raw=True)
        score_cols.append(col)

    cluster_scores = adata.obs.groupby("leiden", observed=True)[score_cols].mean()
    labels = cluster_scores.idxmax(axis=1).str.removeprefix("score_")
    ranked = np.sort(cluster_scores.to_numpy(), axis=1)
    best, second = ranked[:, -1], ranked[:, -2]
    labels[(best < min_score) | (best - second < min_margin)] = UNRESOLVED
    adata.obs["cell_type"] = adata.obs["leiden"].map(labels).astype("category")
    return adata


def find_markers(adata: ad.AnnData, groupby: str, n_genes: int = 50) -> pd.DataFrame:
    """Top genes distinguishing each group from all other cells (Wilcoxon, BH-adjusted)."""
    key = f"rank_genes_{groupby}"
    sc.tl.rank_genes_groups(adata, groupby=groupby, method="wilcoxon", use_raw=True, key_added=key)
    df = sc.get.rank_genes_groups_df(adata, group=None, key=key)
    return df.groupby("group", observed=True).head(n_genes).reset_index(drop=True)


def cd8_mask(adata: ad.AnnData, genes: tuple[str, ...] = ("CD8A", "CD8B")) -> pd.Series:
    """CD8 T cells: CD8A or CD8B detected, outside the clusters that cannot be T cells."""
    present = [g for g in genes if g in adata.raw.var_names]
    if not present:
        raise ValueError(f"none of {genes} are present in this annotation")
    expr = adata.raw[:, present].X
    expr = expr.toarray() if sparse.issparse(expr) else np.asarray(expr)
    detected = (expr > 0).any(axis=1)
    in_t_cluster = ~adata.obs["cell_type"].isin(NON_T_CELL_TYPES).to_numpy()
    return pd.Series(detected & in_t_cluster, index=adata.obs_names)
