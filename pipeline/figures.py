
"""All publication figures. Self-contained: no plotting state outside this module."""
from __future__ import annotations

import pathlib

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.lines import Line2D

# Okabe-Ito first (colour-vision-deficiency safe), then extensions for the long tail.
QUALITATIVE = [
    "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442",
    "#8C6BB1", "#B15928", "#7FC97F", "#386CB0", "#BF5B17", "#666666", "#1B9E77",
    "#D95F02", "#A6761D", "#66A61E", "#E7298A", "#7570B3", "#999999",
]
RESPONSE_COLORS = {"Responder": "#0072B2", "Non-responder": "#D55E00"}
TIMEPOINT_COLORS = {"Pre": "#BFD3E6", "Post": "#2B5D8A"}
GREY = "#4D4D4D"
BASE, SMALL, TINY = 9, 8, 7


def apply_style() -> None:
    mpl.rcParams.update({
        "figure.dpi": 150, "savefig.dpi": 300, "savefig.bbox": "tight",
        "font.size": BASE, "axes.titlesize": BASE, "axes.labelsize": BASE,
        "xtick.labelsize": TINY, "ytick.labelsize": TINY, "legend.fontsize": SMALL,
        "axes.titlelocation": "left", "axes.titleweight": "normal",
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.edgecolor": GREY, "axes.labelcolor": "black",
        "xtick.color": GREY, "ytick.color": GREY,
        "legend.frameon": False, "figure.facecolor": "white",
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def palette_for(categories) -> dict[str, str]:
    cats = [str(c) for c in categories]
    return {c: QUALITATIVE[i % len(QUALITATIVE)] for i, c in enumerate(cats)}


def _embedding_axes(ax, label: str = "UMAP") -> None:
    """§6.6: no ticks on an embedding; a corner arrow pair names the axes."""
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.annotate("", xy=(0.14, 0.02), xytext=(0.02, 0.02), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color=GREY, lw=0.8))
    ax.annotate("", xy=(0.02, 0.14), xytext=(0.02, 0.02), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color=GREY, lw=0.8))
    ax.text(0.155, 0.015, f"{label} 1", transform=ax.transAxes, fontsize=TINY, color=GREY)
    ax.text(0.015, 0.155, f"{label} 2", transform=ax.transAxes, fontsize=TINY,
            color=GREY, rotation=90)


def _scatter(ax, xy, values, colors: dict[str, str], size=2.0, order=None, alpha=0.75):
    values = pd.Series(np.asarray(values).astype(str))
    groups = order or sorted(values.unique())
    for g in groups:
        m = (values == g).to_numpy()
        ax.scatter(xy[m, 0], xy[m, 1], s=size, c=colors.get(g, "#CCCCCC"),
                   linewidths=0, alpha=alpha, rasterized=True)
    ax.margins(0.04)


def _label_on_data(ax, xy, values, fontsize=TINY):
    values = pd.Series(np.asarray(values).astype(str))
    for g in values.unique():
        m = (values == g).to_numpy()
        if m.sum() < 15:
            continue
        x, y = np.median(xy[m, 0]), np.median(xy[m, 1])
        t = ax.text(x, y, g, fontsize=fontsize, ha="center", va="center", color="black")
        t.set_path_effects([mpl.patheffects.withStroke(linewidth=2.2, foreground="white")])


# --------------------------------------------------------------------------- QC
def qc_metrics(obs: pd.DataFrame, qc_summary: pd.DataFrame, cfg, out: pathlib.Path):
    q = cfg.qc
    fig, axes = plt.subplots(1, 4, figsize=(11.5, 2.9))

    panels = [
        ("n_genes", "Genes detected per cell", (q["min_genes"], q["max_genes"]), True),
        ("pct_mito", "Mitochondrial transcripts (%)", (None, q["max_pct_mito"]), False),
        ("total_tpm", "Total TPM assigned per cell", (None, None), True),
    ]
    for ax, (col, label, bounds, logx) in zip(axes, panels):
        vals = obs[col].to_numpy()
        vals = vals[vals > 0] if logx else vals
        bins = np.logspace(np.log10(max(vals.min(), 1)), np.log10(vals.max()), 60) if logx \
            else np.linspace(0, max(vals.max(), 1), 60)
        ax.hist(vals, bins=bins, color="#9EC5E0", edgecolor="white", linewidth=0.2)
        if logx:
            ax.set_xscale("log")
        for b in [x for x in bounds if x is not None]:
            ax.axvline(b, color="#D55E00", lw=1.0, ls="--")
        ax.set_xlabel(label); ax.set_ylabel("Cells")
        ax.margins(0.04)

    ax = axes[3]
    steps = qc_summary["step"].tolist()
    vals = qc_summary["n_cells"].tolist()
    ypos = np.arange(len(steps))[::-1]
    ax.hlines(ypos, 0, vals, color="#B8B8B8", lw=1.0)
    ax.plot(vals, ypos, "o", color="#0072B2", ms=5)
    for v, y in zip(vals, ypos):
        ax.text(v + max(vals) * 0.02, y, f"{v:,}", va="center", fontsize=TINY)
    ax.set_yticks(ypos); ax.set_yticklabels(steps, fontsize=TINY)
    ax.set_xlabel("Cells passing"); ax.set_xlim(0, max(vals) * 1.22)
    ax.set_title("Cells retained at each filter")
    ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)

    n_in = int(qc_summary["n_cells"].iloc[0])
    n_out = int(qc_summary["n_cells"].iloc[-1])
    removed = n_in - n_out
    takeaway = (f"All {n_in:,} deposited cells already satisfy these thresholds "
                f"(dashed) — GEO distributes a pre-filtered matrix"
                if removed == 0 else
                f"{removed:,} of {n_in:,} cells ({100*removed/n_in:.1f}%) fail at least "
                f"one threshold (dashed)")
    fig.suptitle(takeaway, fontsize=BASE, x=0.005, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(out)
    plt.close(fig)
    return out


# ------------------------------------------------------------------------ UMAPs
def umap_clusters(adata, out: pathlib.Path):
    xy = adata.obsm["X_umap"]
    clusters = adata.obs["leiden"].astype(str)
    cpal = palette_for(sorted(clusters.unique(), key=lambda s: int(s)))
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.2))

    _scatter(axes[0], xy, clusters, cpal, order=sorted(cpal, key=int))
    _label_on_data(axes[0], xy, clusters)
    axes[0].set_title(f"Leiden clustering resolves {clusters.nunique()} transcriptional states")

    tp = adata.obs["timepoint"].astype(str)
    shared_tp = int((pd.crosstab(clusters, tp) > 0).all(axis=1).sum())
    _scatter(axes[1], xy, tp, TIMEPOINT_COLORS, order=["Pre", "Post"])
    axes[1].set_title(f"Both timepoints are represented in {shared_tp} of "
                      f"{clusters.nunique()} clusters")
    axes[1].legend(handles=[Line2D([], [], marker="o", ls="", ms=5, color=TIMEPOINT_COLORS[k],
                                   label={"Pre": "Pre-treatment", "Post": "On-treatment"}[k])
                            for k in ["Pre", "Post"]], loc="lower right")

    resp = adata.obs["response"].astype(str)
    shared_r = int((pd.crosstab(clusters, resp) > 0).all(axis=1).sum())
    _scatter(axes[2], xy, resp, RESPONSE_COLORS, order=["Non-responder", "Responder"])
    axes[2].set_title(f"Responders and non-responders both occupy {shared_r} of "
                      f"{clusters.nunique()} clusters")
    axes[2].legend(handles=[Line2D([], [], marker="o", ls="", ms=5, color=v, label=k)
                            for k, v in RESPONSE_COLORS.items()], loc="lower right")

    for ax in axes:
        _embedding_axes(ax)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    return out


def umap_celltypes(adata, out: pathlib.Path):
    xy = adata.obsm["X_umap"]
    ct = adata.obs["cell_type"].astype(str)
    ctpal = palette_for(sorted(ct.unique()))
    ctpal["Unresolved"] = "#D9D9D9"

    state = adata.obs["cd8_state"].astype(str)
    cd8_levels = [s for s in sorted(state.unique()) if s != "Not CD8"]
    spal = {s: QUALITATIVE[i] for i, s in enumerate(cd8_levels)}
    spal["Not CD8"] = "#E8E8E8"

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.0))
    _scatter(axes[0], xy, ct, ctpal, order=sorted(ctpal))
    _label_on_data(axes[0], xy, ct)
    n_lin = len([c for c in ct.unique() if c != "Unresolved"])
    axes[0].set_title(f"{n_lin} canonical lineages recovered from sorted CD45+ cells")

    _scatter(axes[1], xy, state, spal, order=["Not CD8"] + cd8_levels,
             alpha=0.9)
    axes[1].set_title(f"CD8 T cells resolve into {len(cd8_levels)} differentiation states")
    axes[1].legend(handles=[Line2D([], [], marker="o", ls="", ms=5, color=spal[s],
                                   label=s.replace("CD8 ", ""))
                            for s in cd8_levels],
                   loc="center left", bbox_to_anchor=(1.0, 0.5), title="CD8 state")
    for ax in axes:
        _embedding_axes(ax)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    return out


# ---------------------------------------------------------------------- Markers
def marker_dotplot(adata, canonical: dict[str, list[str]], out: pathlib.Path,
                   groupby: str = "population"):
    present = {k: [g for g in v if g in adata.var_names] for k, v in canonical.items()}
    present = {k: v for k, v in present.items() if v}
    sub = adata[adata.obs[groupby].astype(str) != "Unresolved"].copy()
    sub.obs[groupby] = sub.obs[groupby].astype("category").cat.remove_unused_categories()
    dp = sc.pl.dotplot(sub, present, groupby=groupby, standard_scale="var",
                       colorbar_title="Mean expression\n(scaled per gene)",
                       size_title="Fraction of cells (%)", return_fig=True,
                       figsize=(13.0, 5.0))
    dp.style(cmap="Blues", dot_edge_color="white", dot_edge_lw=0.2)
    dp.savefig(out)
    plt.close("all")
    return out


def marker_heatmap(adata, top: dict[str, list[str]], out: pathlib.Path,
                   groupby: str = "population"):
    genes, blocks = [], []
    for pop, gs in top.items():
        genes.extend(gs); blocks.append((pop, len(gs)))
    genes = [g for g in genes if g in adata.var_names]

    sub = adata[adata.obs[groupby].astype(str) != "Unresolved"]
    X = sub[:, genes].X
    X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
    df = pd.DataFrame(X, columns=genes, index=sub.obs[groupby].astype(str).to_numpy())
    mean = df.groupby(level=0).mean()
    z = ((mean - mean.mean(axis=0)) / mean.std(axis=0, ddof=0)).fillna(0.0)
    order = [p for p, _ in blocks if p in z.index]
    z = z.loc[order]

    fig, ax = plt.subplots(figsize=(max(9.0, 0.17 * len(genes)), 0.38 * len(z) + 2.2))
    vmax = float(np.nanpercentile(np.abs(z.to_numpy()), 99))
    im = ax.imshow(z.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=90, fontsize=5.5, style="italic")
    ax.set_yticks(range(len(z))); ax.set_yticklabels(z.index, fontsize=TINY)
    pos = 0
    for pop, n in blocks:
        if pop not in z.index:
            continue
        pos += n
        ax.axvline(pos - 0.5, color="white", lw=1.2)
    ax.set_title(f"Top {max(len(g) for g in top.values())} one-vs-rest markers per "
                 f"population (Wilcoxon, BH q < 0.05)")
    cb = fig.colorbar(im, ax=ax, fraction=0.016, pad=0.012)
    cb.set_label("Mean log$_2$(TPM/10+1), z-scored across populations", fontsize=TINY)
    cb.ax.tick_params(labelsize=TINY)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    return out


# ------------------------------------------------------------------ Composition
def composition_barplots(comp: pd.DataFrame, pop_col: str, out: pathlib.Path):
    pops = sorted(comp[pop_col].astype(str).unique())
    pal = palette_for(pops)
    strata = [("Non-responder", "Pre"), ("Non-responder", "Post"),
              ("Responder", "Pre"), ("Responder", "Post")]

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 4.6),
                             gridspec_kw={"width_ratios": [1, 2.3]})

    ax = axes[0]
    bottoms = np.zeros(len(strata))
    for pop in pops:
        vals = []
        for resp, tp in strata:
            sel = comp[(comp["response"] == resp) & (comp["timepoint"] == tp)
                       & (comp[pop_col].astype(str) == pop)]
            vals.append(sel["fraction"].mean() if len(sel) else 0.0)
        vals = np.asarray(vals)
        ax.bar(range(len(strata)), vals, bottom=bottoms, color=pal[pop], width=0.72,
               edgecolor="white", linewidth=0.3)
        bottoms += vals
    ax.set_xticks(range(len(strata)))
    ax.set_xticklabels([f"{r.replace('Non-responder','NR').replace('Responder','R')}\n{t}"
                        for r, t in strata], fontsize=TINY)
    ax.set_ylabel("Mean fraction of cells per biopsy")
    ax.set_ylim(0, 1); ax.margins(x=0.04)
    ax.set_title("Mean of biopsy compositions, by response and timepoint")

    ax = axes[1]
    order = comp.drop_duplicates("sample").sort_values(["response", "timepoint", "sample"])
    samples = order["sample"].tolist()
    wide = (comp.pivot_table(index="sample", columns=pop_col, values="fraction",
                             observed=True, aggfunc="mean")
            .reindex(samples).fillna(0.0))
    bottoms = np.zeros(len(samples))
    for pop in pops:
        vals = wide[pop].to_numpy() if pop in wide else np.zeros(len(samples))
        ax.bar(range(len(samples)), vals, bottom=bottoms, color=pal[pop], width=0.86,
               edgecolor="white", linewidth=0.15, label=pop)
        bottoms += vals
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=90, fontsize=4.6)
    ax.set_ylim(0, 1); ax.margins(x=0.005)
    ax.set_ylabel("Fraction of cells")
    ax.set_title("")
    ax.legend(loc="center left", bbox_to_anchor=(1.005, 0.5), ncol=1, fontsize=TINY,
              handlelength=0.9, handleheight=0.9, labelspacing=0.28)

    grp = order.groupby(["response", "timepoint"], observed=True).size()
    start = 0
    for (resp, tp), n in grp.items():
        ax.plot([start - 0.4, start + n - 0.6], [1.02, 1.02], color=GREY, lw=1.0,
                clip_on=False)
        ax.text(start + n / 2 - 0.5, 1.04,
                f"{'R' if resp == 'Responder' else 'NR'} {tp}", ha="center",
                fontsize=TINY, color=GREY, clip_on=False)
        start += n
    fig.suptitle("Composition is read per biopsy, then averaged within patient for testing; "
                 "each bar is one biopsy", fontsize=SMALL, x=0.005, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out)
    plt.close(fig)
    return out


def composition_boxplots(pf: pd.DataFrame, stats_df: pd.DataFrame, pop_col: str,
                         out: pathlib.Path, ncol: int = 5):
    pops = stats_df.sort_values("pval")["population"].tolist()
    pops = [p for p in pops if p in set(pf[pop_col].astype(str))]
    nrow = int(np.ceil(len(pops) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.45 * ncol, 2.65 * nrow),
                             squeeze=False)
    rng = np.random.default_rng(0)
    qlook = stats_df.set_index("population")

    for i, pop in enumerate(pops):
        ax = axes[i // ncol][i % ncol]
        sub = pf[pf[pop_col].astype(str) == pop]
        for j, grp in enumerate(["Non-responder", "Responder"]):
            v = sub.loc[sub["response"] == grp, "fraction"].to_numpy()
            if not len(v):
                continue
            bp = ax.boxplot(v, positions=[j], widths=0.5, showfliers=False,
                            patch_artist=True, medianprops=dict(color="black", lw=1.2))
            bp["boxes"][0].set(facecolor=RESPONSE_COLORS[grp], alpha=0.25,
                               edgecolor=RESPONSE_COLORS[grp], linewidth=0.9)
            for part in ("whiskers", "caps"):
                for a in bp[part]:
                    a.set(color=RESPONSE_COLORS[grp], lw=0.9)
            ax.scatter(j + rng.uniform(-0.13, 0.13, len(v)), v, s=11,
                       color=RESPONSE_COLORS[grp], alpha=0.9, linewidths=0, zorder=3)
        row = qlook.loc[pop]
        ax.set_title(f"{pop}\n$p$ = {row['pval']:.3g}   $q$ = {row['qval']:.2g}",
                     fontsize=TINY, loc="center")
        ax.set_xticks([0, 1]); ax.set_xticklabels(["NR", "R"], fontsize=TINY)
        ax.set_ylim(bottom=0); ax.margins(y=0.12)
        if i % ncol == 0:
            ax.set_ylabel("Fraction of patient's cells", fontsize=TINY)
    for k in range(len(pops), nrow * ncol):
        axes[k // ncol][k % ncol].axis("off")
    n_r = pf.drop_duplicates("patient")["response"].value_counts()
    fig.suptitle("Population abundance in responders (R) vs non-responders (NR), "
                 f"one point per patient (n = {int(n_r.get('Responder', 0))} R, "
                 f"{int(n_r.get('Non-responder', 0))} NR); "
                 "Mann-Whitney U, BH-corrected across populations",
                 fontsize=SMALL, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    fig.savefig(out)
    plt.close(fig)
    return out


# ------------------------------------------------------------------- Signatures
def signature_auc(curves: list[dict], out: pathlib.Path):
    fig, ax = plt.subplots(figsize=(5.2, 5.0))
    ax.plot([0, 1], [0, 1], color="#BDBDBD", lw=1.0, ls="--", zorder=1)
    best = max(range(len(curves)), key=lambda j: curves[j]["auc"])
    for i, c in enumerate(curves):
        c = {**c, "focal": i == best}
        ax.plot(c["fpr"], c["tpr"], lw=1.9 if c.get("focal") else 1.3,
                color=QUALITATIVE[i % len(QUALITATIVE)],
                alpha=1.0 if c.get("focal") else 0.85,
                label=f"{c['name']} — AUC {c['auc']:.2f} (p = {c['pval']:.3f})", zorder=3)
    ax.set_xlabel("False-positive rate (non-responders called responder)")
    ax.set_ylabel("True-positive rate (responders recovered)")
    ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)
    ax.set_aspect("equal")
    ax.set_title("Leave-one-patient-out discrimination of response")
    ax.legend(loc="lower right", fontsize=TINY)
    ax.text(0.03, 0.97, "chance = diagonal", transform=ax.transAxes, fontsize=TINY,
            color=GREY, va="top")
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    return out


def signature_scores(patient_scores: pd.DataFrame, cell_scores: pd.DataFrame,
                     out: pathlib.Path):
    panels = [c for c in patient_scores.columns if c != "response"]
    fig, axes = plt.subplots(1, len(panels) + 1, figsize=(2.9 * (len(panels) + 1), 3.4),
                             squeeze=False)
    axes = axes[0]
    rng = np.random.default_rng(0)
    for ax, col in zip(axes, panels):
        for j, grp in enumerate(["Non-responder", "Responder"]):
            v = patient_scores.loc[patient_scores["response"] == grp, col].dropna().to_numpy()
            ax.scatter(j + rng.uniform(-0.12, 0.12, len(v)), v, s=22,
                       color=RESPONSE_COLORS[grp], alpha=0.9, linewidths=0)
            ax.hlines(np.median(v), j - 0.26, j + 0.26, color="black", lw=1.6, zorder=4)
        ax.set_xticks([0, 1]); ax.set_xticklabels(["NR", "R"])
        ax.set_xlim(-0.6, 1.6); ax.margins(y=0.14)
        ax.set_title(col, fontsize=SMALL)
        ax.set_ylabel("Patient score" if ax is axes[0] else "")

    ax = axes[-1]
    for grp in ["Non-responder", "Responder"]:
        v = cell_scores.loc[cell_scores["response"] == grp, "score"].to_numpy()
        ax.hist(v, bins=60, density=True, histtype="stepfilled", alpha=0.42,
                color=RESPONSE_COLORS[grp], label=grp)
        ax.hist(v, bins=60, density=True, histtype="step", lw=1.1,
                color=RESPONSE_COLORS[grp])
    ax.set_xlabel("Per-cell signature score"); ax.set_ylabel("Density")
    ax.set_title("Single CD8 T cells", fontsize=SMALL)
    ax.legend(loc="upper right", fontsize=TINY)
    fig.suptitle("Signature scores separate responders at the patient level; "
                 "single cells overlap heavily", fontsize=SMALL, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out)
    plt.close(fig)
    return out
