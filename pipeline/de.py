"""Patient-level pseudobulk differential expression, responder vs non-responder.

Cells from one patient are not independent replicates, so a cell-level test
would report p-values driven by cell counts rather than by patients. Every
contrast here is computed on patient pseudobulk: the mean log2(TPM/10+1) of a
gene over that patient's cells in the population of interest.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats
from statsmodels.stats.multitest import multipletests


def pseudobulk(adata, cfg, population: str | None = None, by: str = "patient"):
    """Return (expression genes x units, unit metadata, detection rate per gene).

    Aggregation is on the TPM layer -- the mean TPM over a patient's cells -- and
    the log transform is applied afterwards. Averaging log values instead would
    compress every fold-change towards zero, because most cells contribute a log
    of a zero-inflated value rather than a zero-inflated value itself.
    """
    sub = adata if population is None else adata[adata.obs["population"] == population]
    d = cfg["differential_expression"]
    scale = float(cfg["normalisation"]["tpm_scale"])

    X = sub.layers["tpm"] if "tpm" in sub.layers else sub.X
    X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)
    units = sub.obs[by].astype(str).to_numpy()
    keep_units = [u for u, n in pd.Series(units).value_counts().items()
                  if n >= d["min_cells_per_patient"]]
    keep_units = sorted(keep_units)

    means, ncells = {}, {}
    for u in keep_units:
        idx = np.flatnonzero(units == u)
        means[u] = np.asarray(X[idx].mean(axis=0)).ravel()
        ncells[u] = len(idx)
    expr = pd.DataFrame(means, index=sub.var_names)[keep_units]
    expr = np.log2(expr / scale + 1.0)       # pseudobulk TPM -> log2(TPM/10+1)

    detect = np.asarray((X > 0).mean(axis=0)).ravel()
    detection = pd.Series(detect, index=sub.var_names, name="detection_rate")

    meta = (sub.obs.groupby(by, observed=True)[["response", "therapy"]]
            .agg(lambda s: s.iloc[0]).loc[keep_units])
    meta["n_cells"] = pd.Series(ncells)
    meta["responder"] = (meta["response"] == "Responder").astype(int)
    return expr, meta, detection


def filter_genes(expr: pd.DataFrame, detection: pd.Series, cfg) -> pd.DataFrame:
    d = cfg["differential_expression"]
    keep = (detection >= d["min_detection_rate"]) & \
           ((expr > 0).sum(axis=1) >= d["min_patients_expressing"])
    return expr.loc[keep.reindex(expr.index, fill_value=False)]


def test_response(expr: pd.DataFrame, meta: pd.DataFrame, cfg,
                  population: str = "all") -> pd.DataFrame:
    """Welch t-test per gene on patient pseudobulk; BH-FDR across genes."""
    y = meta["responder"].to_numpy()
    R = expr.loc[:, meta.index[y == 1]].to_numpy()
    N = expr.loc[:, meta.index[y == 0]].to_numpy()
    t, p = stats.ttest_ind(R, N, axis=1, equal_var=False)

    auc = np.array([_auc(expr.iloc[i].to_numpy(), y) for i in range(expr.shape[0])])
    out = pd.DataFrame({
        "population": population,
        "gene": expr.index,
        "mean_log2_R": R.mean(axis=1),
        "mean_log2_NR": N.mean(axis=1),
        "log2FC": R.mean(axis=1) - N.mean(axis=1),   # inputs are already log2
        "t_stat": t,
        "auc": auc,
        "pval": p,
        "n_R": int((y == 1).sum()),
        "n_NR": int((y == 0).sum()),
    })
    out["qval"] = multipletests(out["pval"].fillna(1.0), method="fdr_bh")[1]
    out["significant"] = out["qval"] < cfg["differential_expression"]["fdr_alpha"]
    return out.sort_values("pval").reset_index(drop=True)


def _auc(x: np.ndarray, y: np.ndarray) -> float:
    pos, neg = x[y == 1], x[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    u = stats.mannwhitneyu(pos, neg, alternative="two-sided").statistic
    return float(u / (len(pos) * len(neg)))
