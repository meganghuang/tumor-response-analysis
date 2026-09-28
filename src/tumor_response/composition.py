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


CONFLICTING = "Conflicting"


def _single_label(values: pd.Series) -> str:
    """The one label shared by a patient's samples, or `Conflicting` if they disagree."""
    unique = values.dropna().unique()
    return str(unique[0]) if len(unique) == 1 else CONFLICTING


def aggregate_to_patient(
    comp: pd.DataFrame,
    cell_types: list[str],
    patient_col: str = "patient",
    annotations: tuple[str, ...] = ("response", "therapy"),
) -> pd.DataFrame:
    """One row per patient: mean cell-type fraction across that patient's samples.

    Samples are not independent -- 13 of the 32 patients contribute more than one biopsy --
    so the patient is the correct unit for a between-patient comparison. Response is
    annotated per biopsy in this dataset and four patients have lesions that disagree; those
    patients are labelled `Conflicting` rather than being silently assigned to a group.
    """
    out = comp.groupby(patient_col)[cell_types].mean()
    for col in annotations:
        if col in comp.columns:
            out[col] = comp.groupby(patient_col)[col].agg(_single_label)
    out.index.name = patient_col
    out[patient_col] = out.index.astype(str)
    return out


def conflicting_patients(
    comp: pd.DataFrame, patient_col: str = "patient", response_col: str = "response"
) -> list[str]:
    """Patients whose biopsies carry more than one distinct response label."""
    n = comp.groupby(patient_col)[response_col].nunique()
    return sorted(str(p) for p in n[n > 1].index)


def therapy_switchers(
    comp: pd.DataFrame, patient_col: str = "patient", therapy_col: str = "therapy"
) -> list[str]:
    """Patients whose biopsies were taken under more than one therapy."""
    n = comp.groupby(patient_col)[therapy_col].nunique()
    return sorted(str(p) for p in n[n > 1].index)


def single_timepoint_samples(
    comp: pd.DataFrame, patient_col: str = "patient", time_col: str = "timepoint"
) -> pd.DataFrame:
    """Samples from patients biopsied at only one timepoint.

    An unpaired Pre-versus-Post test needs two independent groups, but 11 patients here
    contribute to both. Restricting to patients present at a single timepoint makes the two
    groups genuinely disjoint; the paired test covers the rest.
    """
    n_time = comp.groupby(patient_col)[time_col].nunique()
    keep = n_time[n_time == 1].index
    return comp[comp[patient_col].isin(keep)]


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
