import anndata as ad
import scanpy as sc


def annotate_qc_metrics(adata: ad.AnnData) -> ad.AnnData:
    """Add per-cell QC metrics, including percent mitochondrial reads."""
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    return adata


def filter_cells_and_genes(
    adata: ad.AnnData,
    min_genes: int,
    max_genes: int,
    min_cells: int,
    max_pct_mito: float,
) -> ad.AnnData:
    """Remove low-quality cells, likely doublets, and rarely detected genes."""
    if "pct_counts_mt" not in adata.obs:
        annotate_qc_metrics(adata)
    keep = (
        (adata.obs["n_genes_by_counts"] >= min_genes)
        & (adata.obs["n_genes_by_counts"] <= max_genes)
        & (adata.obs["pct_counts_mt"] <= max_pct_mito)
    )
    adata = adata[keep].copy()
    sc.pp.filter_genes(adata, min_cells=min_cells)
    return adata
