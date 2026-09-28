"""Quantitative comparison of the scanpy and Seurat paths.

Both paths read the same raw GEO files and use the same QC thresholds, marker sets and
statistics, so the remaining differences are implementation choices: variable-gene selection
(scanpy's dispersion-based `seurat` flavor vs Seurat's `vst`), community detection (Leiden vs
Louvain), UMAP initialisation, module scoring (`score_genes` vs `AddModuleScore`), and the
marker-gene pre-filter. This script measures how much those choices move the results.
"""
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score

from tumor_response.config import load_config, project_path
from tumor_response.plotting import save, set_style


def jaccard(a: set, b: set) -> float:
    return len(a & b) / len(a | b) if (a | b) else np.nan


def main():
    cfg = load_config()
    processed_dir = project_path(cfg["data"]["processed_dir"])
    tables_dir = project_path(cfg["output"]["tables_dir"])
    figures_dir = project_path(cfg["output"]["figures_dir"])
    seurat_tables = project_path(f"{cfg['output']['tables_dir']}/seurat")
    set_style()

    adata = sc.read_h5ad(processed_dir / "02_annotated.h5ad")
    py = pd.DataFrame({"cluster": adata.obs["leiden"].astype(str),
                       "cell_type": adata.obs["cell_type"].astype(str),
                       "sample": adata.obs[cfg["data"]["sample_column"]].astype(str)},
                      index=adata.obs_names)
    r = pd.read_csv(seurat_tables / "cell_labels.csv", index_col="cell")

    shared = py.index.intersection(r.index)
    rows = [{"metric": "cells in scanpy path", "value": len(py)},
            {"metric": "cells in Seurat path", "value": len(r)},
            {"metric": "cells in both", "value": len(shared)}]

    # Clustering agreement on the cells both paths kept.
    a, b = py.loc[shared, "cluster"], r.loc[shared, "cluster"].astype(str)
    rows += [
        {"metric": "clusters, scanpy (Leiden)", "value": a.nunique()},
        {"metric": "clusters, Seurat (Louvain)", "value": b.nunique()},
        {"metric": "adjusted Rand index, clusters", "value": adjusted_rand_score(a, b)},
        {"metric": "adjusted mutual information, clusters",
         "value": adjusted_mutual_info_score(a, b)},
    ]

    # Cell-type label agreement.
    ta, tb = py.loc[shared, "cell_type"], r.loc[shared, "cell_type"]
    rows += [
        {"metric": "cell types, scanpy", "value": ta.nunique()},
        {"metric": "cell types, Seurat", "value": tb.nunique()},
        {"metric": "cell-type labels identical (fraction of cells)", "value": (ta == tb).mean()},
        {"metric": "adjusted Rand index, cell types", "value": adjusted_rand_score(ta, tb)},
    ]
    pd.crosstab(ta, tb).to_csv(tables_dir / "concordance_cell_type_confusion.csv")

    # Per-sample composition correlation, per cell type shared by both paths.
    comp_py = pd.read_csv(tables_dir / "composition_by_sample.csv", index_col=0)
    comp_r = pd.read_csv(seurat_tables / "composition_by_sample.csv", index_col="sample")
    # Only the cell-type columns; both tables also carry sample annotations.
    cell_types = set(adata.obs["cell_type"].astype(str)) | set(r["cell_type"].astype(str))
    shared_types = [c for c in comp_py.columns if c in comp_r.columns and c in cell_types]
    comp_rows = []
    for ct in shared_types:
        x = comp_py.loc[comp_py.index.intersection(comp_r.index), ct]
        y = comp_r.loc[x.index, ct]
        comp_rows.append({"cell_type": ct, "n_samples": len(x),
                          "pearson_r": stats.pearsonr(x, y)[0] if x.std() and y.std() else np.nan,
                          "spearman_r": stats.spearmanr(x, y)[0] if x.std() and y.std() else np.nan,
                          "mean_fraction_scanpy": x.mean(), "mean_fraction_seurat": y.mean()})
    comp_table = pd.DataFrame(comp_rows).set_index("cell_type")
    comp_table.to_csv(tables_dir / "concordance_composition.csv")
    rows.append({"metric": "median per-sample composition Pearson r (shared cell types)",
                 "value": comp_table["pearson_r"].median()})

    # Marker-gene overlap per cell type, top 50 each side. Only available when the Seurat
    # path was run with TUMOR_RESPONSE_MARKERS=1; see R/README.md for why it is off by
    # default (no presto on this platform, so Seurat's Wilcoxon is prohibitively slow).
    jac_table = None
    r_markers = seurat_tables / "markers_by_cell_type.csv"
    if r_markers.exists():
        mk_py = pd.read_csv(tables_dir / "markers_by_cell_type.csv")
        mk_r = pd.read_csv(r_markers)
        jac = []
        for ct in sorted(set(mk_py["group"].astype(str)) & set(mk_r["cluster"].astype(str))):
            top_py = set(mk_py.loc[mk_py["group"].astype(str) == ct, "names"].head(50))
            top_r = set(mk_r.loc[mk_r["cluster"].astype(str) == ct, "gene"].head(50))
            jac.append({"cell_type": ct, "n_scanpy": len(top_py), "n_seurat": len(top_r),
                        "shared": len(top_py & top_r), "jaccard": jaccard(top_py, top_r)})
        jac_table = pd.DataFrame(jac).set_index("cell_type")
        jac_table.to_csv(tables_dir / "concordance_markers.csv")
        rows.append({"metric": "median top-50 marker Jaccard (shared cell types)",
                     "value": jac_table["jaccard"].median()})
    else:
        rows.append({"metric": "marker overlap", "value": np.nan})

    summary = pd.DataFrame(rows).set_index("metric")
    summary.to_csv(tables_dir / "concordance_summary.csv")

    # Headline AUCs side by side.
    auc_py = pd.read_csv(tables_dir / "response_signature_summary.csv")
    auc_r = pd.read_csv(seurat_tables / "response_signature_summary.csv")
    key = ["signature", "unit"]
    auc = auc_py.merge(auc_r, on=key, how="outer", suffixes=("_scanpy", "_seurat"))
    auc = auc[key + ["n_scanpy", "auc_scanpy", "n_seurat", "auc_seurat"]]
    auc["auc_difference"] = auc["auc_seurat"] - auc["auc_scanpy"]
    auc.to_csv(tables_dir / "concordance_auc.csv", index=False)

    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
    ax = axes[0]
    ax.scatter(comp_table["mean_fraction_scanpy"], comp_table["mean_fraction_seurat"],
               s=30, color="#2a78d6")
    lim = max(comp_table[["mean_fraction_scanpy", "mean_fraction_seurat"]].max()) * 1.1
    ax.plot([0, lim], [0, lim], color="#6b6a64", lw=0.8, ls="--")
    for ct, row in comp_table.iterrows():
        ax.annotate(ct, (row["mean_fraction_scanpy"], row["mean_fraction_seurat"]), fontsize=6,
                    xytext=(3, 2), textcoords="offset points")
    ax.set_xlabel("mean fraction, scanpy")
    ax.set_ylabel("mean fraction, Seurat")
    ax.set_title("Cell-type abundance", loc="left")

    ax = axes[1]
    if jac_table is not None:
        ok = jac_table.dropna(subset=["jaccard"]).sort_values("jaccard")
        ax.barh(range(len(ok)), ok["jaccard"], color="#eb6834")
        ax.set_yticks(range(len(ok)), ok.index, fontsize=7)
        ax.set_xlabel("Jaccard overlap, top 50 markers")
        ax.set_title("Marker agreement", loc="left")
    else:
        conf = pd.crosstab(ta, tb, normalize="index")
        im = ax.imshow(conf.to_numpy(), cmap="Blues", vmin=0, vmax=1, aspect="auto")
        ax.set_xticks(range(conf.shape[1]), conf.columns, rotation=40, ha="right", fontsize=6)
        ax.set_yticks(range(conf.shape[0]), conf.index, fontsize=6)
        ax.set_xlabel("Seurat label")
        ax.set_ylabel("scanpy label")
        ax.set_title("Cell-type label agreement", loc="left")
        ax.grid(False)
        fig.colorbar(im, ax=ax, fraction=0.04, label="fraction of scanpy label")

    ax = axes[2]
    valid = auc.dropna(subset=["auc_scanpy", "auc_seurat"])
    ax.scatter(valid["auc_scanpy"], valid["auc_seurat"], s=40, color="#2a78d6")
    ax.plot([0.3, 1.0], [0.3, 1.0], color="#6b6a64", lw=0.8, ls="--")
    ax.axhline(0.5, color="#e4e3dd", lw=0.8)
    ax.axvline(0.5, color="#e4e3dd", lw=0.8)
    ax.set_xlabel("AUC, scanpy")
    ax.set_ylabel("AUC, Seurat")
    ax.set_title("Response signatures", loc="left")
    fig.suptitle("scanpy vs Seurat, same data and same statistics", x=0.01, ha="left")
    fig.tight_layout()
    save(fig, figures_dir / "concordance_python_vs_seurat.png")

    print(summary.round(3).to_string())
    print(auc.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
