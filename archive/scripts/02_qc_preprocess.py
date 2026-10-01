"""Load GSE120575, run QC filtering, and compute PCA / neighbors / UMAP."""
import scanpy as sc

from tumor_response.config import load_config, project_path
from tumor_response.io import load_sade_feldman
from tumor_response.preprocess import normalize_and_reduce
from tumor_response.qc import annotate_qc_metrics, filter_cells_and_genes


def main():
    cfg = load_config()
    raw_dir = project_path(cfg["data"]["raw_dir"])
    processed_dir = project_path(cfg["data"]["processed_dir"])
    sc.settings.figdir = project_path(cfg["output"]["figures_dir"])

    # Parsing the 55k-gene text matrix takes several minutes, so cache it once.
    cache = processed_dir / "00_raw.h5ad"
    if cache.exists():
        adata = sc.read_h5ad(cache)
    else:
        adata = load_sade_feldman(raw_dir)
        adata.write_h5ad(cache)
    print(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes, {adata.obs['sample'].nunique()} samples")

    log_normalized = cfg["data"]["already_log_normalized"]
    annotate_qc_metrics(adata, log_normalized=log_normalized)
    print(f"Mitochondrial genes in this annotation: {adata.uns['n_mito_genes']}")
    metrics = ["n_genes_by_counts", "pct_counts_mt"]
    if "pct_counts_mt_logscale" in adata.obs:
        metrics.append("pct_counts_mt_logscale")
        print("pct mito, linear scale vs uncorrected log scale (mean, 95th pct): "
              f"{adata.obs['pct_counts_mt'].mean():.2f} / "
              f"{adata.obs['pct_counts_mt'].quantile(0.95):.2f} vs "
              f"{adata.obs['pct_counts_mt_logscale'].mean():.2f} / "
              f"{adata.obs['pct_counts_mt_logscale'].quantile(0.95):.2f}")
    sc.pl.violin(adata, metrics, multi_panel=True, save="_qc.png", show=False)

    q = cfg["qc"]
    adata = filter_cells_and_genes(adata, q["min_genes"], q["max_genes"], q["min_cells"],
                                   q["max_pct_mito"], log_normalized=log_normalized)
    print(f"Cells removed by each QC rule: {adata.uns['qc_removed']}")
    print(f"After QC: {adata.n_obs} cells x {adata.n_vars} genes")

    p = cfg["preprocess"]
    adata = normalize_and_reduce(adata, p["target_sum"], p["n_top_genes"], p["n_pcs"],
                                 cfg["clustering"]["n_neighbors"],
                                 already_log_normalized=cfg["data"]["already_log_normalized"])
    adata.write_h5ad(processed_dir / "01_preprocessed.h5ad")


if __name__ == "__main__":
    main()
