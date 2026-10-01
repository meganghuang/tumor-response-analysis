"""Unit tests on synthetic data: statistics and annotation logic, not biology."""
from __future__ import annotations

import sys
import pathlib

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from pipeline import annotate, composition, config, de, embed, qc, signatures


@pytest.fixture(scope="module")
def cfg():
    return config.load_config()


def _toy(n_cells=400, n_genes=60, seed=0):
    rng = np.random.default_rng(seed)
    X = rng.poisson(5, size=(n_cells, n_genes)).astype(np.float32) * 20.0
    genes = [f"G{i}" for i in range(n_genes - 4)] + ["MT-A", "MT-B", "CD8A", "CD8B"]
    patients = rng.choice([f"P{i}" for i in range(8)], n_cells)
    obs = pd.DataFrame({
        "sample": [f"{'Pre' if i % 2 else 'Post'}_{p}" for i, p in enumerate(patients)],
        "patient": patients,
        "timepoint": ["Pre" if i % 2 else "Post" for i in range(n_cells)],
        "response": np.where(pd.Series(patients).isin([f"P{i}" for i in range(4)]),
                             "Responder", "Non-responder"),
        "therapy": "anti-PD1",
    }, index=[f"c{i}" for i in range(n_cells)])
    a = ad.AnnData(sp.csr_matrix(X), obs=obs, var=pd.DataFrame(index=genes))
    for c in ["sample", "patient", "timepoint", "response", "therapy"]:
        a.obs[c] = a.obs[c].astype("category")
    return a


def test_qc_filters_and_reports_each_step(cfg):
    a = _toy()
    a.X[0, :] = 0.0                      # an empty well
    a, summary = qc.apply_qc(a, cfg)
    assert "n_genes" in a.obs and "pct_mito" in a.obs
    assert summary["n_cells"].is_monotonic_decreasing or summary["n_cells"].iloc[-1] <=         summary["n_cells"].iloc[0]
    assert a.n_obs < 400


def test_normalise_is_log2_of_tpm_over_scale(cfg):
    a = _toy()
    raw = a.X.copy()
    embed.normalise(a, cfg)
    expected = np.log2(raw[:5].toarray() / cfg["normalisation"]["tpm_scale"] + 1)
    np.testing.assert_allclose(a.X[:5].toarray(), expected, rtol=1e-5)
    assert "tpm" in a.layers


def test_assign_labels_marks_unresolved_only_when_nothing_is_above_average():
    z = pd.DataFrame({"A": [2.0, -1.0], "B": [-1.0, -2.0]}, index=["c0", "c1"])
    labels = annotate.assign_labels(z, min_z=0.0)
    assert labels["c0"] == "A"
    assert labels["c1"] == "Unresolved"


def test_paired_treatment_test_recovers_a_planted_shift(cfg):
    rows = []
    for i in range(12):
        rows += [{"patient": f"P{i}", "timepoint": "Pre", "population": "X",
                  "fraction": 0.2, "sample": f"Pre_P{i}", "n_cells": 100,
                  "response": "Responder", "therapy": "a", "n_total": 500},
                 {"patient": f"P{i}", "timepoint": "Post", "population": "X",
                  "fraction": 0.4, "sample": f"Post_P{i}", "n_cells": 200,
                  "response": "Responder", "therapy": "a", "n_total": 500}]
    out = composition.treatment_paired_test(pd.DataFrame(rows), "population", 0.1)
    assert out.loc[0, "n_increase"] == 12
    assert out.loc[0, "log2FC"] > 0.9
    assert out.loc[0, "pval"] < 0.01


def test_response_test_auc_matches_perfect_separation():
    rows = [{"patient": f"P{i}", "population": "X", "fraction": 0.1 + 0.01 * i,
             "n_cells": 50, "n_samples": 1,
             "response": "Responder" if i >= 10 else "Non-responder"} for i in range(20)]
    out = composition.response_test(pd.DataFrame(rows), "population", 0.1, "resp")
    assert out.loc[0, "auc"] == 1.0
    assert out.loc[0, "pval"] < 0.001


def test_pseudobulk_collapses_cells_to_patients(cfg):
    a = _toy()
    embed.normalise(a, cfg)
    expr, meta, detection = de.pseudobulk(a, cfg)
    assert expr.shape[1] == len(meta) <= 8
    assert set(expr.columns) == set(meta.index)
    assert (detection <= 1).all()


def test_lopo_never_uses_the_held_out_patient():
    rng = np.random.default_rng(0)
    genes = [f"G{i}" for i in range(50)]
    patients = [f"P{i}" for i in range(10)]
    expr = pd.DataFrame(rng.normal(size=(50, 10)), index=genes, columns=patients)
    meta = pd.DataFrame({"responder": [1] * 5 + [0] * 5}, index=patients)
    expr.iloc[0, :5] += 5.0                      # one strongly discriminative gene
    s = signatures.lopo_scores(expr, meta, size=3)
    assert len(s) == 10 and s.notna().all()
    # Shuffling a held-out patient's own label must not change its own score.
    meta2 = meta.copy(); meta2.loc["P0", "responder"] = 0
    s2 = signatures.lopo_scores(expr, meta2, size=3)
    assert np.isclose(s["P0"], s2["P0"])


def test_permutation_p_is_bounded_and_auc_is_one_for_perfect_scores():
    y = np.array([0] * 8 + [1] * 8)
    score = np.arange(16, dtype=float)
    auc, p = signatures.auc_with_permutation(y, score, n_perm=200, seed=0)
    assert auc == 1.0
    assert 0 < p <= 1
