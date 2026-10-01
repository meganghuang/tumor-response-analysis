"""Marker-programme annotation of lineages and CD8 T-cell states.

Clusters are labelled by scoring every cell against a panel of canonical
programmes, averaging within cluster, z-scoring each programme across clusters
and taking the argmax. A cluster is Unresolved only if no programme is above
average anywhere -- deliberately avoiding a winner-take-all margin rule, which
discards the T-cell clusters where two lineage programmes are legitimately close.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc

LINEAGE_MARKERS: dict[str, list[str]] = {
    "CD8 T": ["CD8A", "CD8B", "CD3D", "CD3E", "CD3G", "GZMK"],
    "CD4 T": ["CD4", "IL7R", "CD40LG", "CD3D", "CD3E", "TNFRSF4", "CCR7"],
    "Treg": ["FOXP3", "IL2RA", "IKZF2", "CTLA4", "TNFRSF9", "TNFRSF18"],
    "NK": ["KLRF1", "NCAM1", "FGFBP2", "FCGR3A", "KLRD1", "GNLY", "NKG7", "TYROBP"],
    "B": ["MS4A1", "CD79A", "CD79B", "CD19", "BANK1", "TCL1A"],
    "Plasma": ["MZB1", "DERL3", "XBP1", "IGHG1", "IGKC", "TNFRSF17"],
    "Mono/Macrophage": ["LYZ", "CD14", "CD68", "CSF1R", "AIF1", "C1QA", "C1QB", "APOE"],
    "DC": ["CD1C", "FCER1A", "CLEC9A", "CLEC10A", "CST3", "HLA-DQA1"],
    "pDC": ["LILRA4", "CLEC4C", "IL3RA", "IRF7", "TCF4", "SERPINF1"],
    "Cycling": ["MKI67", "TOP2A", "CCNB1", "UBE2C", "TYMS", "STMN1", "BIRC5"],
    "Melanocytic": ["MLANA", "PMEL", "TYR", "DCT", "MITF", "SOX10"],
}

CD8_STATE_MARKERS: dict[str, list[str]] = {
    "CD8 memory-like (TCF7+)": ["TCF7", "LEF1", "SELL", "CCR7", "IL7R", "BACH2", "FOXP1"],
    "CD8 effector/cytotoxic": ["GZMB", "GZMH", "PRF1", "FGFBP2", "KLRG1", "CX3CR1", "NKG7"],
    "CD8 exhausted/dysfunctional": ["PDCD1", "LAG3", "HAVCR2", "TIGIT", "CTLA4", "ENTPD1",
                                    "TOX", "CD38", "LAYN"],
    "CD8 tissue-resident": ["ITGAE", "CD69", "ZNF683", "XCL1", "XCL2", "RGS1"],
    "CD8 cycling": ["MKI67", "TOP2A", "UBE2C", "BIRC5", "TYMS"],
}


def score_programmes(adata, programmes: dict[str, list[str]], prefix: str = "score_") -> list[str]:
    """Score every cell against each programme with ``sc.tl.score_genes``."""
    keys = []
    for name, genes in programmes.items():
        present = [g for g in genes if g in adata.var_names]
        if len(present) < 2:
            continue
        key = prefix + name
        sc.tl.score_genes(adata, present, score_name=key, random_state=0, use_raw=False)
        keys.append(key)
    return keys


def cluster_programme_matrix(adata, keys: list[str], cluster_key: str) -> pd.DataFrame:
    """Mean programme score per cluster, z-scored across clusters within programme."""
    df = adata.obs.groupby(cluster_key, observed=True)[keys].mean()
    z = (df - df.mean(axis=0)) / df.std(axis=0, ddof=0).replace(0, np.nan)
    z.columns = [c.split("score_", 1)[1] for c in z.columns]
    return z


def assign_labels(zmat: pd.DataFrame, min_z: float) -> pd.Series:
    best = zmat.idxmax(axis=1)
    best[zmat.max(axis=1) <= min_z] = "Unresolved"
    return best


def annotate_lineages(adata, cfg, cluster_key: str = "leiden"):
    """Label Leiden clusters with lineage identities. Returns the z-score matrix."""
    keys = score_programmes(adata, LINEAGE_MARKERS)
    zmat = cluster_programme_matrix(adata, keys, cluster_key)
    labels = assign_labels(zmat, cfg["annotation"]["min_z"])
    adata.obs["cell_type"] = adata.obs[cluster_key].map(labels).astype("category")
    return zmat


def annotate_cd8_states(adata, sub, cfg, sub_key: str = "cd8_leiden"):
    """Label CD8 subclusters with differentiation states and write back to ``adata``."""
    keys = score_programmes(sub, CD8_STATE_MARKERS)
    zmat = cluster_programme_matrix(sub, keys, sub_key)
    labels = assign_labels(zmat, cfg["annotation"]["min_z"])
    sub.obs["cd8_state"] = sub.obs[sub_key].map(labels)

    state = pd.Series("Not CD8", index=adata.obs_names, dtype=object)
    state.loc[sub.obs_names] = sub.obs["cd8_state"].values
    adata.obs["cd8_state"] = pd.Categorical(state)

    pop = adata.obs["cell_type"].astype(str).copy()
    pop.loc[sub.obs_names] = sub.obs["cd8_state"].astype(str).values
    adata.obs["population"] = pd.Categorical(pop)
    return zmat
