"""Population composition: expansion/contraction with treatment and with response.

Two questions, two units of analysis:

* **Treatment effect** -- paired within patient. Patients contributing both a
  pre- and an on-treatment biopsy are compared with a Wilcoxon signed-rank test
  on their population fractions.
* **Response association** -- the patient is the unit, never the cell and never
  the biopsy, because biopsies from one patient are not independent. Fractions
  are averaged across a patient's biopsies within the relevant timepoint stratum
  and compared between responders and non-responders with Mann-Whitney U.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


def sample_composition(adata, pop_col: str = "population", unit: str = "sample",
                       min_cells: int = 0) -> pd.DataFrame:
    """Long table of population fractions per biopsy, with clinical annotation."""
    obs = adata.obs
    counts = pd.crosstab(obs[unit], obs[pop_col])
    counts = counts.loc[counts.sum(axis=1) >= min_cells]
    frac = counts.div(counts.sum(axis=1), axis=0)
    long = (frac.stack().rename("fraction").reset_index()
            .merge(counts.stack().rename("n_cells").reset_index(),
                   on=[unit, pop_col]))
    meta = (obs.groupby(unit, observed=True)[["patient", "timepoint", "response", "therapy"]]
            .agg(lambda s: s.iloc[0]).reset_index())
    long = long.merge(meta, on=unit, how="left")
    long["n_total"] = long[unit].map(counts.sum(axis=1))
    return long


def patient_fractions(comp: pd.DataFrame, pop_col: str, timepoint: str | None = None,
                      unit: str = "sample") -> pd.DataFrame:
    """Collapse biopsies to one value per patient (mean over that patient's biopsies)."""
    df = comp if timepoint is None else comp[comp["timepoint"] == timepoint]
    out = (df.groupby(["patient", pop_col], observed=True)
             .agg(fraction=("fraction", "mean"), n_cells=("n_cells", "sum"),
                  n_samples=(unit, "nunique"))
             .reset_index())
    resp = df.groupby("patient", observed=True)["response"].agg(lambda s: s.iloc[0])
    out["response"] = out["patient"].map(resp)
    return out


def _auc(x: np.ndarray, y: np.ndarray) -> float:
    """AUC of ``x`` separating group 1 from group 0 in ``y`` (Mann-Whitney identity)."""
    pos, neg = x[y == 1], x[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    u = stats.mannwhitneyu(pos, neg, alternative="two-sided").statistic
    return float(u / (len(pos) * len(neg)))


def treatment_paired_test(comp: pd.DataFrame, pop_col: str, alpha: float) -> pd.DataFrame:
    """Wilcoxon signed-rank on pre vs on-treatment fractions, paired within patient."""
    wide = (comp.groupby(["patient", "timepoint", pop_col], observed=True)["fraction"]
            .mean().unstack("timepoint"))
    if not {"Pre", "Post"}.issubset(wide.columns):
        raise ValueError("both Pre and Post timepoints are required")
    wide = wide.dropna(subset=["Pre", "Post"]).reset_index()

    rows = []
    for pop, grp in wide.groupby(pop_col, observed=True):
        pre, post = grp["Pre"].to_numpy(), grp["Post"].to_numpy()
        delta = post - pre
        try:
            pval = 1.0 if np.allclose(delta, 0) else float(
                stats.wilcoxon(post, pre, zero_method="wilcox").pvalue)
        except ValueError:
            pval = 1.0
        rows.append({
            "analysis": "treatment_paired", "population": str(pop), "n_patients": len(grp),
            "mean_Pre": pre.mean(), "mean_Post": post.mean(),
            "median_delta": float(np.median(delta)),
            "log2FC": float(np.log2((post.mean() + 1e-4) / (pre.mean() + 1e-4))),
            "n_increase": int((delta > 0).sum()), "n_decrease": int((delta < 0).sum()),
            "statistic": np.nan, "pval": pval,
        })
    return _fdr(pd.DataFrame(rows), alpha)


def response_test(pf: pd.DataFrame, pop_col: str, alpha: float, label: str) -> pd.DataFrame:
    """Mann-Whitney U on patient-level fractions, responder vs non-responder."""
    rows = []
    for pop, grp in pf.groupby(pop_col, observed=True):
        y = (grp["response"] == "Responder").astype(int).to_numpy()
        x = grp["fraction"].to_numpy()
        if y.sum() < 3 or (1 - y).sum() < 3:
            continue
        res = stats.mannwhitneyu(x[y == 1], x[y == 0], alternative="two-sided")
        rows.append({
            "analysis": label, "population": str(pop),
            "n_responder": int(y.sum()), "n_nonresponder": int((1 - y).sum()),
            "mean_R": float(x[y == 1].mean()), "mean_NR": float(x[y == 0].mean()),
            "median_R": float(np.median(x[y == 1])), "median_NR": float(np.median(x[y == 0])),
            "log2FC": float(np.log2((x[y == 1].mean() + 1e-4) / (x[y == 0].mean() + 1e-4))),
            "auc": _auc(x, y), "statistic": float(res.statistic), "pval": float(res.pvalue),
        })
    return _fdr(pd.DataFrame(rows), alpha)


def _fdr(df: pd.DataFrame, alpha: float) -> pd.DataFrame:
    if df.empty:
        return df
    df = df.copy()
    df["qval"] = multipletests(df["pval"], method="fdr_bh")[1]
    df["significant"] = df["qval"] < alpha
    return df.sort_values("pval").reset_index(drop=True)
