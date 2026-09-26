import anndata as ad
import numpy as np
import pandas as pd

from tumor_response.qc import annotate_qc_metrics, filter_cells_and_genes


def make_adata():
    rng = np.random.default_rng(0)
    X = rng.poisson(1.0, size=(50, 30)).astype(np.float32)
    X[0, :] = 0  # empty cell
    X[1, 0] = 500  # cell dominated by a mitochondrial gene
    genes = ["MT-CO1"] + [f"GENE{i}" for i in range(1, 30)]
    return ad.AnnData(X, var=pd.DataFrame(index=genes))


def test_annotate_qc_metrics_flags_mito_genes():
    adata = annotate_qc_metrics(make_adata())
    assert adata.var.loc["MT-CO1", "mt"]
    assert adata.var["mt"].sum() == 1
    assert "pct_counts_mt" in adata.obs


def test_filter_removes_empty_and_high_mito_cells():
    adata = make_adata()
    filtered = filter_cells_and_genes(adata, min_genes=5, max_genes=1000, min_cells=1, max_pct_mito=50)
    kept = set(filtered.obs_names)
    assert "0" not in kept
    assert "1" not in kept
    assert filtered.n_obs == 48
