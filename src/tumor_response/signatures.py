"""Response signatures: sample-level (pseudobulk) differential expression and scoring.

Cells from the same sample are not independent, so group comparisons are done on one value
per sample (the sample's mean expression) rather than on individual cells.
"""
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats

from tumor_response.composition import _bh


def pseudobulk_matrices(
    adata: ad.AnnData,
    sample_col: str,
    mask: pd.Series | None = None,
    min_cells: int = 10,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.Series]:
    """Per-sample mean log expression, detection rate, and cell count, over all genes.

    Only cells in `mask` are used, and samples with fewer than `min_cells` such cells are
    dropped. No gene filtering is applied here: the detection rates are returned so that a
    detection filter can be recomputed from any subset of samples -- which is what
    cross-validation needs, since a filter fitted on all samples has seen the held-out
    patient.
    """
    raw = adata.raw.to_adata()
    if mask is not None:
        raw = raw[mask.loc[raw.obs_names].to_numpy()]
    X = sparse.csr_matrix(raw.X)

    samples = raw.obs[sample_col].astype(str)
    sizes = samples.value_counts()
    kept = sizes.index[sizes >= min_cells]
    onehot = sparse.csr_matrix(pd.get_dummies(samples)[kept].to_numpy(dtype=np.float32))
    n_cells = pd.Series(np.asarray(onehot.sum(axis=0)).ravel(), index=kept, name="n_cells")

    sums = np.asarray((onehot.T @ X).todense())
    detected = np.asarray((onehot.T @ (X > 0).astype(np.float32)).todense())
    denom = n_cells.to_numpy()[:, None]
    means = pd.DataFrame(sums / denom, index=kept, columns=raw.var_names)
    detect = pd.DataFrame(detected / denom, index=kept, columns=raw.var_names)
    return means, detect, n_cells


def detection_rate(
    detect: pd.DataFrame, n_cells: pd.Series, samples: pd.Index | None = None
) -> pd.Series:
    """Fraction of cells with the gene detected, pooled over `samples` (default: all).

    Weighting each sample by its cell count reproduces the pooled per-cell detection rate.
    """
    if samples is not None:
        detect, n_cells = detect.loc[samples], n_cells.loc[samples]
    return detect.mul(n_cells, axis=0).sum() / n_cells.sum()


def pseudobulk(
    adata: ad.AnnData,
    sample_col: str,
    mask: pd.Series | None = None,
    min_cells: int = 10,
    min_frac_expressed: float = 0.1,
) -> pd.DataFrame:
    """Mean log expression (from `adata.raw`) per sample, samples x genes.

    Genes detected in fewer than `min_frac_expressed` of the cells are dropped. Use
    `pseudobulk_matrices` instead when the gene filter must be refitted per fold.
    """
    means, detect, n_cells = pseudobulk_matrices(adata, sample_col, mask, min_cells)
    keep = detection_rate(detect, n_cells) >= min_frac_expressed
    return means.loc[:, keep.to_numpy()]


def pseudobulk_de(pb: pd.DataFrame, labels: pd.Series, group_a: str, group_b: str) -> pd.DataFrame:
    """Per-gene Mann-Whitney U test between samples of two groups.

    `log2fc` is the difference in mean ln(TPM+1) converted to log2 units (group_b vs group_a);
    `auc` is how well the gene alone separates the groups (>0.5: higher in group_b).
    """
    labels = labels.loc[pb.index]
    a, b = pb[labels == group_a], pb[labels == group_b]
    test = stats.mannwhitneyu(b.to_numpy(), a.to_numpy(), axis=0, alternative="two-sided")
    out = pd.DataFrame(
        {
            f"mean_{group_a}": a.mean(),
            f"mean_{group_b}": b.mean(),
            "log2fc": (b.mean() - a.mean()) / np.log(2),
            "auc": test.statistic / (len(a) * len(b)),
            "pvalue": test.pvalue,
        },
        index=pb.columns,
    )
    out["qvalue"] = _bh(out["pvalue"])
    out.index.name = "gene"
    return out.sort_values("pvalue")


def top_genes(de: pd.DataFrame, n: int) -> tuple[list[str], list[str]]:
    """The n most significant genes up and down in group_b.

    Fewer than `n` are returned if fewer genes move in that direction, so callers that
    report a signature size should report the lengths actually obtained.
    """
    up = de[de["log2fc"] > 0].head(n).index.tolist()
    down = de[de["log2fc"] < 0].head(n).index.tolist()
    return up, down


def by_patient(values: pd.Series, patients: pd.Series) -> pd.Series:
    """Average a per-sample value within each patient, giving one value per patient."""
    return values.groupby(patients.loc[values.index].to_numpy()).mean()


def signature_score(pb: pd.DataFrame, up: list[str], down: list[str],
                    ref: pd.DataFrame | None = None) -> pd.Series:
    """Mean z-score of `up` genes minus mean z-score of `down` genes, per sample.

    z-scores use the mean/std of `ref` (default: `pb` itself), so held-out samples can be
    scored with parameters learned only from training samples.
    """
    ref = pb if ref is None else ref
    z = (pb - ref.mean()) / ref.std().replace(0, np.nan)
    score = z[up].mean(axis=1) if up else 0
    if down:
        score = score - z[down].mean(axis=1)
    return score


def auc(scores: pd.Series, labels: pd.Series, positive: str) -> float:
    """ROC AUC: P(score of a random `positive` sample > score of a random other sample)."""
    labels = labels.loc[scores.index]
    pos, neg = scores[labels == positive], scores[labels != positive]
    return stats.mannwhitneyu(pos, neg).statistic / (len(pos) * len(neg))


def leave_one_patient_out_scores(
    pb: pd.DataFrame, labels: pd.Series, patients: pd.Series,
    group_a: str, group_b: str, n_genes: int,
    detect: pd.DataFrame | None = None,
    n_cells: pd.Series | None = None,
    min_frac_expressed: float = 0.0,
) -> pd.Series:
    """Score every sample with a signature learned without any sample from its patient.

    For each patient, the signature (top `n_genes` up/down genes) is derived from the other
    patients' samples and applied to the held-out patient's samples. The AUC of these scores
    is an honest estimate of how well such a signature generalizes to new patients.

    Everything fitted is fitted on training samples only: the detection filter (when
    `detect`/`n_cells` are supplied), the differential expression that picks the genes, and
    the mean/std used to z-score. Pass `pb` unfiltered so the per-fold filter is meaningful.
    """
    labels, patients = labels.loc[pb.index], patients.loc[pb.index]
    scores = pd.Series(np.nan, index=pb.index)
    for patient in patients.unique():
        test = patients == patient
        train = ~test
        if labels[train].nunique() < 2:
            continue
        sub = pb
        if detect is not None and n_cells is not None and min_frac_expressed > 0:
            keep = detection_rate(detect, n_cells, pb.index[train]) >= min_frac_expressed
            sub = pb.loc[:, keep.reindex(pb.columns).fillna(False).to_numpy()]
        de = pseudobulk_de(sub[train], labels[train], group_a, group_b)
        up, down = top_genes(de, n_genes)
        scores[test] = signature_score(sub[test], up, down, ref=sub[train])
    return scores


def fraction_positive(adata: ad.AnnData, gene: str, mask: pd.Series, sample_col: str,
                      min_cells: int = 10) -> pd.Series:
    """Per sample, the fraction of cells in `mask` with the gene detected (e.g. TCF7+ of CD8)."""
    cells = mask[mask].index
    x = adata.raw[cells, [gene]].X
    x = x.toarray() if sparse.issparse(x) else np.asarray(x)
    df = pd.DataFrame({"sample": adata.obs.loc[cells, sample_col].astype(str).values,
                       "pos": x.ravel() > 0})
    grouped = df.groupby("sample")["pos"]
    frac = grouped.mean()
    return frac[grouped.size() >= min_cells]
