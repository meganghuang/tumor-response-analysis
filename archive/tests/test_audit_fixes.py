"""Tests for the corrections listed in AUDIT.md."""
import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from tumor_response.composition import (
    CONFLICTING, aggregate_to_patient, composition_by_sample, conflicting_patients,
    single_timepoint_samples, therapy_switchers,
)
from tumor_response.qc import annotate_qc_metrics
from tumor_response.signatures import (
    detection_rate, leave_one_patient_out_scores, pseudobulk_matrices, pseudobulk_de, top_genes,
)


# --- A1: mitochondrial fraction on the linear scale ---------------------------------------

def make_log_adata():
    """One cell dominated by a mitochondrial gene, stored as ln(TPM+1)."""
    counts = np.array([
        [900.0, 50.0, 50.0],   # 90% mitochondrial
        [10.0, 495.0, 495.0],  # 1% mitochondrial
    ], dtype=np.float32)
    adata = ad.AnnData(np.log1p(counts), var=pd.DataFrame(index=["MT-CO1", "GENE1", "GENE2"]))
    adata.uns["true_pct"] = 100 * counts[:, 0] / counts.sum(axis=1)
    return adata


def test_mito_percent_recovers_the_linear_scale_fraction():
    adata = annotate_qc_metrics(make_log_adata(), log_normalized=True)
    assert np.allclose(adata.obs["pct_counts_mt"], adata.uns["true_pct"], atol=1e-3)


def test_uncorrected_log_scale_fraction_understates_high_mito_cells():
    adata = annotate_qc_metrics(make_log_adata(), log_normalized=True)
    high = adata.obs.iloc[0]
    assert high["pct_counts_mt"] == pytest.approx(90.0, abs=1e-3)
    # The log-scale ratio the original code reported is far lower, which is why a 20%
    # threshold let high-mitochondrial cells through.
    assert high["pct_counts_mt_logscale"] < 50
    assert adata.uns["n_mito_genes"] == 1


def test_mito_percent_is_zero_without_mito_genes():
    adata = ad.AnnData(np.log1p(np.ones((3, 2), dtype=np.float32)),
                       var=pd.DataFrame(index=["GENE1", "GENE2"]))
    adata = annotate_qc_metrics(adata, log_normalized=True)
    assert adata.uns["n_mito_genes"] == 0
    assert (adata.obs["pct_counts_mt"] == 0).all()


def test_sparse_and_dense_agree():
    dense = annotate_qc_metrics(make_log_adata(), log_normalized=True)
    sp = make_log_adata()
    sp.X = sparse.csr_matrix(sp.X)
    sp = annotate_qc_metrics(sp, log_normalized=True)
    assert np.allclose(dense.obs["pct_counts_mt"], sp.obs["pct_counts_mt"])


# --- A2/A3: patient as the unit, and conflicting response labels --------------------------

def make_comp():
    """Four patients. P1 has two biopsies that disagree on response; P2 switches therapy."""
    rows = [
        ("Pre_P1", "P1", "Pre", "Responder", "anti-CTLA4", 0.8),
        ("Post_P1", "P1", "Post", "Non-responder", "anti-PD1", 0.6),
        ("Pre_P2", "P2", "Pre", "Responder", "anti-CTLA4", 0.7),
        ("Post_P2", "P2", "Post", "Responder", "anti-PD1", 0.9),
        ("Pre_P3", "P3", "Pre", "Non-responder", "anti-PD1", 0.2),
        ("Post_P4", "P4", "Post", "Non-responder", "anti-PD1", 0.1),
    ]
    df = pd.DataFrame(rows, columns=["sample", "patient", "timepoint", "response", "therapy",
                                     "memory"]).set_index("sample")
    df["exhausted"] = 1 - df["memory"]
    return df


def test_conflicting_patients_are_detected():
    assert conflicting_patients(make_comp()) == ["P1"]


def test_therapy_switchers_are_detected():
    assert therapy_switchers(make_comp()) == ["P1", "P2"]


def test_aggregate_to_patient_averages_and_marks_conflicts():
    out = aggregate_to_patient(make_comp(), ["memory", "exhausted"])
    assert len(out) == 4
    assert out.loc["P2", "memory"] == pytest.approx(0.8)  # mean of 0.7 and 0.9
    assert out.loc["P1", "response"] == CONFLICTING
    assert out.loc["P2", "response"] == "Responder"
    assert out.loc["P2", "therapy"] == CONFLICTING  # two therapies, one patient


def test_conflicted_patients_drop_out_of_both_groups():
    out = aggregate_to_patient(make_comp(), ["memory", "exhausted"])
    grouped = out[out["response"].isin(["Responder", "Non-responder"])]
    assert set(grouped.index) == {"P2", "P3", "P4"}


# --- A4: independent groups for the unpaired timepoint test -------------------------------

def test_single_timepoint_samples_drops_patients_present_at_both():
    kept = single_timepoint_samples(make_comp())
    assert set(kept.index) == {"Pre_P3", "Post_P4"}
    assert kept.groupby("patient")["timepoint"].nunique().max() == 1


# --- A5: the detection filter is refitted inside each fold --------------------------------

def make_adata(n_per_sample=20, seed=0):
    """Six patients, two biopsies each. GENE0 tracks response; GENE1 exists only in P0."""
    rng = np.random.default_rng(seed)
    obs, rows = [], []
    for p in range(6):
        resp = "Responder" if p < 3 else "Non-responder"
        for tp in ["Pre", "Post"]:
            for _ in range(n_per_sample):
                obs.append({"sample": f"{tp}_P{p}", "patient": f"P{p}", "timepoint": tp,
                            "response": resp, "therapy": "anti-PD1"})
                x = rng.poisson(1.0, 10).astype(np.float32)
                x[0] += 5 if resp == "Responder" else 0
                x[1] = rng.poisson(5.0) if p == 0 else 0.0
                rows.append(np.log1p(x))
    obs = pd.DataFrame(obs, index=[str(i) for i in range(len(obs))]).astype("category")
    adata = ad.AnnData(sparse.csr_matrix(np.array(rows)), obs=obs,
                       var=pd.DataFrame(index=[f"GENE{i}" for i in range(10)]))
    adata.raw = adata
    return adata


def test_detection_rate_uses_only_the_named_samples():
    means, detect, n_cells = pseudobulk_matrices(make_adata(), "sample", min_cells=5)
    assert means.shape == (12, 10)
    overall = detection_rate(detect, n_cells)
    p0_samples = [s for s in means.index if s.endswith("P0")]
    without_p0 = detection_rate(detect, n_cells, means.index.difference(p0_samples))
    # GENE1 is present only in P0, so it clears a global filter but not a filter fitted
    # without that patient -- the leak the fix removes.
    assert overall["GENE1"] > 0.1
    assert without_p0["GENE1"] == 0.0
    # GENE0 is expressed in every patient, so dropping one barely moves its detection rate.
    assert overall["GENE0"] > 0.5 and without_p0["GENE0"] > 0.5
    assert abs(overall["GENE0"] - without_p0["GENE0"]) < 0.1


def test_detection_rate_weights_samples_by_cell_count():
    detect = pd.DataFrame({"G": [1.0, 0.0]}, index=["a", "b"])
    n_cells = pd.Series([30, 10], index=["a", "b"])
    assert detection_rate(detect, n_cells)["G"] == pytest.approx(0.75)


def test_in_fold_filter_keeps_the_planted_signal_and_scores_every_sample():
    adata = make_adata()
    means, detect, n_cells = pseudobulk_matrices(adata, "sample", min_cells=5)
    info = adata.obs.groupby("sample", observed=True)[["response", "patient"]].first().astype(str)
    scores = leave_one_patient_out_scores(
        means, info["response"], info["patient"], "Non-responder", "Responder", n_genes=1,
        detect=detect, n_cells=n_cells, min_frac_expressed=0.1)
    assert scores.notna().all()
    de = pseudobulk_de(means, info["response"], "Non-responder", "Responder")
    assert de.index[0] == "GENE0"


# --- C3: a short signature is reported honestly -------------------------------------------

def test_top_genes_returns_fewer_when_few_genes_move_that_way():
    de = pd.DataFrame({"log2fc": [1.0, 2.0, -0.5], "pvalue": [0.01, 0.02, 0.03]},
                      index=["A", "B", "C"])
    up, down = top_genes(de, n=5)
    assert up == ["A", "B"] and down == ["C"]


def test_composition_by_sample_still_sums_to_one():
    adata = make_adata()
    adata.obs["cell_type"] = pd.Categorical(
        np.where(np.arange(adata.n_obs) % 2 == 0, "memory", "exhausted"))
    comp = composition_by_sample(adata, "sample")
    assert np.allclose(comp[["memory", "exhausted"]].sum(axis=1), 1)
