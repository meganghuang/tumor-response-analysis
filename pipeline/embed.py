"""Normalisation, feature selection, dimensionality reduction and clustering.

Scaling is applied to a highly-variable-gene view only. The full log-expression
matrix is left untouched and sparse, so every downstream test, score and plot
reads the same log2(TPM/10+1) values rather than a z-scored copy.
"""
from __future__ import annotations

import numpy as np
import scanpy as sc
import scipy.sparse as sp


def normalise(adata, cfg) -> None:
    """log2(TPM/scale + 1), following Sade-Feldman et al. Raw TPM kept in a layer."""
    scale = float(cfg["normalisation"]["tpm_scale"])
    adata.layers["tpm"] = adata.X.copy()
    X = adata.X / scale
    if sp.issparse(X):
        X.data = np.log2(X.data + 1.0)
    else:
        X = np.log2(X + 1.0)
    adata.X = X


def _pca_graph_cluster(adata, cfg, cluster_key: str, resolution: float, umap: bool):
    e = cfg.embedding
    sc.pp.highly_variable_genes(adata, n_top_genes=e["n_top_genes"], flavor=e["hvg_flavor"])
    hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(hvg, max_value=e["scale_max_value"], zero_center=True)
    n_comps = int(min(e["n_pcs"], hvg.n_obs - 1, hvg.n_vars - 1))
    sc.tl.pca(hvg, n_comps=n_comps, svd_solver="arpack", random_state=e["random_state"])

    adata.obsm["X_pca"] = hvg.obsm["X_pca"]
    adata.uns["pca"] = hvg.uns["pca"]
    del hvg

    # Exact kNN. The approximate (pynndescent) backend spawns a joblib worker
    # pool, which is unavailable in restricted environments; at this cohort size
    # the exact graph costs seconds and is deterministic.
    sc.pp.neighbors(adata, n_neighbors=e["n_neighbors"], n_pcs=n_comps,
                    random_state=e["random_state"],
                    transformer=e.get("knn_transformer", "sklearn"))
    sc.tl.leiden(adata, resolution=resolution, key_added=cluster_key,
                 random_state=e["random_state"], flavor="igraph", n_iterations=2,
                 directed=False)
    if umap:
        sc.tl.umap(adata, random_state=e["random_state"])
    return adata


def embed(adata, cfg):
    """HVGs -> PCA -> kNN -> Leiden -> UMAP on the whole dataset, in place."""
    return _pca_graph_cluster(adata, cfg, "leiden",
                              cfg.embedding["leiden_resolution"], umap=True)


def subcluster(adata, mask, cfg, key: str, resolution: float):
    """Re-embed a subset of cells and Leiden-cluster it. Returns the subset AnnData."""
    sub = adata[mask].copy()
    return _pca_graph_cluster(sub, cfg, key, resolution, umap=True)
