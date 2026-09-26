import anndata as ad
import scanpy as sc


def normalize_and_reduce(
    adata: ad.AnnData,
    target_sum: float,
    n_top_genes: int,
    n_pcs: int,
    n_neighbors: int,
    already_log_normalized: bool = False,
) -> ad.AnnData:
    """Normalize, select variable genes, and compute PCA, neighbors, and UMAP.

    Set `already_log_normalized` for data that is already ln(x+1)-scaled (e.g. TPM from
    Smart-seq2); otherwise raw counts are kept in `adata.layers["counts"]` and normalized.
    Log-normalized values are always kept in `adata.raw`.
    """
    if not already_log_normalized:
        adata.layers["counts"] = adata.X.copy()
        sc.pp.normalize_total(adata, target_sum=target_sum)
        sc.pp.log1p(adata)
    adata.raw = adata
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=n_pcs, mask_var="highly_variable")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.umap(adata)
    return adata
