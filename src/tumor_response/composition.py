"""Cell-type composition per sample and group comparisons."""
import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats


def composition_by_sample(
    adata: ad.AnnData,
    sample_col: str,
    cell_type_col: str = "cell_type",
    annotations: tuple[str, ...] = ("patient", "timepoint", "response", "therapy"),
) -> pd.DataFrame:
    """Fraction of each cell type per sample, plus sample-level annotations.

    Proportions are computed per sample (not per cell) so samples with more cells
    don't dominate the comparison. Cell-type columns come first; annotations last.
    """
    counts = pd.crosstab(adata.obs[sample_col], adata.obs[cell_type_col])
    fractions = counts.div(counts.sum(axis=1), axis=0)
    fractions.columns = fractions.columns.astype(str)
    sample_info = adata.obs.groupby(sample_col, observed=True)[list(annotations)].first()
    return fractions.join(sample_info.astype(str))


def _bh(pvalues: pd.Series) -> pd.Series:
    q = np.full(len(pvalues), np.nan)
    ok = pvalues.notna().to_numpy()
    if ok.any():
        q[ok] = stats.false_discovery_control(pvalues[ok].to_numpy())
    return pd.Series(q, index=pvalues.index)


def compare_fractions(
    comp: pd.DataFrame, cell_types: list[str], group_col: str, group_a: str, group_b: str
) -> pd.DataFrame:
    """Unpaired comparison (Mann-Whitney U) of each cell type's fraction between two groups.

    `auc` is P(fraction in a random group_b sample > one in group_a); 0.5 means no difference.
    """
    a = comp[comp[group_col] == group_a]
    b = comp[comp[group_col] == group_b]
    rows = []
    for ct in cell_types:
        test = stats.mannwhitneyu(b[ct], a[ct], alternative="two-sided")
        rows.append({
            "cell_type": ct,
            f"mean_{group_a}": a[ct].mean(),
            f"mean_{group_b}": b[ct].mean(),
            "diff": b[ct].mean() - a[ct].mean(),
            "auc": test.statistic / (len(a) * len(b)),
            "pvalue": test.pvalue,
        })
    out = pd.DataFrame(rows).set_index("cell_type")
    out["qvalue"] = _bh(out["pvalue"])
    out.attrs["n"] = {group_a: len(a), group_b: len(b)}
    return out.sort_values("pvalue")


def paired_fractions(
    comp: pd.DataFrame,
    cell_types: list[str],
    patient_col: str = "patient",
    time_col: str = "timepoint",
    before: str = "Pre",
    after: str = "Post",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Paired before/after comparison (Wilcoxon signed-rank) for patients sampled at both times.

    Patients with several lesions at one timepoint are averaged. Returns the test table and
    the per-patient fractions used (long format) for plotting.
    """
    per_patient = comp.groupby([patient_col, time_col])[cell_types].mean()
    wide = per_patient.unstack(time_col)
    paired = wide.dropna(how="any")
    rows = []
    for ct in cell_types:
        pre, post = paired[(ct, before)], paired[(ct, after)]
        diff = post - pre
        test = stats.wilcoxon(post, pre) if (diff != 0).any() else None
        rows.append({
            "cell_type": ct,
            f"mean_{before}": pre.mean(),
            f"mean_{after}": post.mean(),
            "mean_change": diff.mean(),
            "n_increase": int((diff > 0).sum()),
            "n_decrease": int((diff < 0).sum()),
            "pvalue": test.pvalue if test else np.nan,
        })
    out = pd.DataFrame(rows).set_index("cell_type")
    out["qvalue"] = _bh(out["pvalue"])
    out.attrs["n_patients"] = len(paired)
    long = paired.stack(time_col, future_stack=True).reset_index()
    return out.sort_values("pvalue"), long
