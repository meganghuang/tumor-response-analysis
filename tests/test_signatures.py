import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from tumor_response.composition import compare_fractions, composition_by_sample, paired_fractions
from tumor_response.signatures import (
    auc, fraction_positive, leave_one_patient_out_scores, pseudobulk, pseudobulk_de,
)


def make_adata(n_per_sample=30, seed=0):
    """8 patients x (Pre, Post); responders have more "memory" cells and express GENE0 higher."""
    rng = np.random.default_rng(seed)
    obs, rows = [], []
    for p in range(8):
        resp = "Responder" if p < 4 else "Non-responder"
        for tp in ["Pre", "Post"]:
            for _ in range(n_per_sample):
                is_memory = rng.random() < (0.7 if resp == "Responder" else 0.2)
                obs.append({"sample": f"{tp}_P{p}", "patient": f"P{p}", "timepoint": tp,
                            "response": resp, "therapy": "anti-PD1",
                            "cell_type": "memory" if is_memory else "exhausted"})
                x = rng.poisson(1.0, 20).astype(np.float32)
                x[0] += 5 if resp == "Responder" else 0
                rows.append(np.log1p(x))
    obs = pd.DataFrame(obs, index=[str(i) for i in range(len(obs))]).astype("category")
    adata = ad.AnnData(sparse.csr_matrix(np.array(rows)), obs=obs,
                       var=pd.DataFrame(index=[f"GENE{i}" for i in range(20)]))
    adata.raw = adata
    return adata


def test_composition_sums_to_one_and_detects_enrichment():
    adata = make_adata()
    comp = composition_by_sample(adata, "sample")
    assert np.allclose(comp[["memory", "exhausted"]].sum(axis=1), 1)
    table = compare_fractions(comp, ["memory", "exhausted"], "response", "Non-responder", "Responder")
    assert table.loc["memory", "diff"] > 0
    assert table.loc["memory", "auc"] > 0.9


def test_paired_fractions_uses_patients_with_both_timepoints():
    comp = composition_by_sample(make_adata(), "sample")
    comp = comp.drop(index="Post_P0")
    table, long = paired_fractions(comp, ["memory"])
    assert table.attrs["n_patients"] == 7
    assert len(long) == 14


def test_pseudobulk_de_finds_planted_gene_and_signature_generalizes():
    adata = make_adata()
    pb = pseudobulk(adata, "sample", min_cells=5, min_frac_expressed=0)
    assert pb.shape == (16, 20)
    info = adata.obs.groupby("sample", observed=True)[["response", "patient"]].first().astype(str)
    de = pseudobulk_de(pb, info["response"], "Non-responder", "Responder")
    assert de.index[0] == "GENE0" and de.loc["GENE0", "log2fc"] > 0
    scores = leave_one_patient_out_scores(pb, info["response"], info["patient"],
                                          "Non-responder", "Responder", n_genes=1)
    assert auc(scores, info["response"], "Responder") > 0.9


def test_fraction_positive_counts_only_masked_cells():
    adata = make_adata()
    mask = pd.Series(adata.obs["cell_type"] == "memory", index=adata.obs_names)
    frac = fraction_positive(adata, "GENE0", mask, "sample", min_cells=1)
    assert frac.between(0, 1).all()
    assert (frac[frac.index.str.contains("P0")] == 1).all()  # responders always express GENE0
