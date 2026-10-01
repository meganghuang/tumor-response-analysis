"""Per-cell quality control.

The mitochondrial fraction has to be computed on the linear (TPM) scale. This dataset
arrives as log2(TPM+1) and is converted to ln(TPM+1) on load, so summing `X` directly --
which is what `scanpy.pp.calculate_qc_metrics` does -- gives a ratio of sums of logarithms
rather than a fraction of transcripts. Because the logarithm compresses large values, that
ratio is biased downward for exactly the high-mitochondrial cells the filter exists to
remove.
"""
import anndata as ad
import numpy as np
import scanpy as sc
from scipy import sparse


def mito_percent(adata: ad.AnnData) -> np.ndarray:
    """Percent of TPM contributed by mitochondrial genes, computed on the linear scale.

    Assumes `adata.X` holds ln(TPM+1); values are back-transformed with expm1 before the
    fraction is taken. Returns 0 for cells with no detected expression.
    """
    mt = adata.var["mt"].to_numpy()
    X = adata.X
    if sparse.issparse(X):
        linear = X.tocsr(copy=True).astype(np.float64)
        linear.data = np.expm1(linear.data)
    else:
        linear = np.expm1(np.asarray(X, dtype=np.float64))
    total = np.asarray(linear.sum(axis=1)).ravel()
    mito = np.asarray(linear[:, mt].sum(axis=1)).ravel() if mt.any() else np.zeros_like(total)
    return np.where(total > 0, 100.0 * mito / total, 0.0)


def annotate_qc_metrics(adata: ad.AnnData, log_normalized: bool = False) -> ad.AnnData:
    """Add per-cell QC metrics, including percent mitochondrial expression.

    Set `log_normalized` when `adata.X` is already ln(x+1)-scaled, so the mitochondrial
    fraction is computed on back-transformed values. The uncorrected log-scale ratio is kept
    as `pct_counts_mt_logscale` for comparison, and the number of mitochondrial genes found
    is recorded in `adata.uns` so a missing `MT-` annotation cannot pass unnoticed.
    """
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    adata.uns["n_mito_genes"] = int(adata.var["mt"].sum())
    if log_normalized:
        adata.obs["pct_counts_mt_logscale"] = adata.obs["pct_counts_mt"].to_numpy()
        adata.obs["pct_counts_mt"] = mito_percent(adata)
    return adata


def filter_cells_and_genes(
    adata: ad.AnnData,
    min_genes: int,
    max_genes: int,
    min_cells: int,
    max_pct_mito: float,
    log_normalized: bool = False,
) -> ad.AnnData:
    """Remove low-quality cells and rarely detected genes.

    `max_genes` is kept for continuity with droplet-data conventions, but these are
    plate-based Smart-seq2 libraries with one cell per well, so a high detected-gene count is
    not doublet evidence; how many cells it removes is recorded in `adata.uns["qc_removed"]`.
    """
    if "pct_counts_mt" not in adata.obs:
        annotate_qc_metrics(adata, log_normalized=log_normalized)
    n_genes = adata.obs["n_genes_by_counts"]
    reasons = {
        "too_few_genes": int((n_genes < min_genes).sum()),
        "too_many_genes": int((n_genes > max_genes).sum()),
        "high_mito": int((adata.obs["pct_counts_mt"] > max_pct_mito).sum()),
    }
    keep = (
        (n_genes >= min_genes)
        & (n_genes <= max_genes)
        & (adata.obs["pct_counts_mt"] <= max_pct_mito)
    )
    removed = reasons | {"total_cells_removed": int((~keep).sum())}
    adata = adata[keep].copy()
    sc.pp.filter_genes(adata, min_cells=min_cells)
    adata.uns["qc_removed"] = removed
    return adata
