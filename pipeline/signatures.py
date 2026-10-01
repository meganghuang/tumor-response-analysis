"""Learn and evaluate transcriptional signatures that stratify responders.

A signature derived from all patients and then scored on those same patients is
circular. Every headline AUC reported here is leave-one-patient-out: the gene
list and the gene-wise standardisation are re-derived from the training patients
only, and the held-out patient is scored with them. The all-patient gene list is
still exported, as the object a future cohort would be scored with.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
from sklearn.metrics import roc_auc_score, roc_curve


def rank_by_t(expr: pd.DataFrame, y: np.ndarray) -> pd.Series:
    R = expr.loc[:, y == 1].to_numpy()
    N = expr.loc[:, y == 0].to_numpy()
    t, _ = stats.ttest_ind(R, N, axis=1, equal_var=False)
    return pd.Series(np.nan_to_num(t), index=expr.index).sort_values(ascending=False)


def build_signature(expr: pd.DataFrame, y: np.ndarray, size: int) -> tuple[list[str], list[str]]:
    """Top ``size`` genes up in responders and up in non-responders."""
    t = rank_by_t(expr, y)
    return list(t.index[:size]), list(t.index[-size:][::-1])


def score_units(expr: pd.DataFrame, up: list[str], down: list[str],
                centre: pd.Series, scale: pd.Series) -> pd.Series:
    """Mean standardised expression of up genes minus that of down genes."""
    z = expr.sub(centre, axis=0).div(scale.replace(0, np.nan), axis=0)
    return z.loc[up].mean(axis=0) - z.loc[down].mean(axis=0)


def lopo_scores(expr: pd.DataFrame, meta: pd.DataFrame, size: int) -> pd.Series:
    """Leave-one-patient-out signature score for every patient."""
    y_all = meta["responder"].to_numpy()
    scores = {}
    for i, patient in enumerate(meta.index):
        train = np.ones(len(meta), bool); train[i] = False
        tr_expr = expr.loc[:, meta.index[train]]
        up, down = build_signature(tr_expr, y_all[train], size)
        centre, scale = tr_expr.mean(axis=1), tr_expr.std(axis=1, ddof=1)
        scores[patient] = float(
            score_units(expr.loc[:, [patient]], up, down, centre, scale).iloc[0]
        )
    return pd.Series(scores, name="lopo_score")


def auc_with_permutation(y: np.ndarray, score: np.ndarray, n_perm: int,
                         seed: int) -> tuple[float, float]:
    """Observed AUC and a two-sided permutation p-value against shuffled labels."""
    ok = np.isfinite(score)
    y, score = y[ok], score[ok]
    obs = roc_auc_score(y, score)
    rng = np.random.default_rng(seed)
    null = np.array([roc_auc_score(rng.permutation(y), score) for _ in range(n_perm)])
    p = (np.sum(np.abs(null - 0.5) >= abs(obs - 0.5)) + 1) / (n_perm + 1)
    return float(obs), float(p)


def roc_points(y: np.ndarray, score: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    ok = np.isfinite(score)
    fpr, tpr, _ = roc_curve(y[ok], score[ok])
    return fpr, tpr


def score_cells(adata, up: list[str], down: list[str], name: str) -> str:
    """Per-cell signature score (up programme minus down programme)."""
    up = [g for g in up if g in adata.var_names]
    down = [g for g in down if g in adata.var_names]
    sc.tl.score_genes(adata, up, score_name=f"_{name}_up", random_state=0, use_raw=False)
    sc.tl.score_genes(adata, down, score_name=f"_{name}_dn", random_state=0, use_raw=False)
    adata.obs[name] = adata.obs[f"_{name}_up"] - adata.obs[f"_{name}_dn"]
    del adata.obs[f"_{name}_up"], adata.obs[f"_{name}_dn"]
    return name
