"""Shared figure style and plots for two-group comparisons."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Colorblind-safe categorical pair (validated for CVD separation); used for every
# two-group comparison: first group, second group.
PAIR = ["#2a78d6", "#eb6834"]
INK, MUTED, GRID = "#1f1f1e", "#6b6a64", "#e4e3dd"


def set_style():
    plt.rcParams.update({
        "figure.dpi": 110, "savefig.dpi": 200, "savefig.bbox": "tight",
        "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
        "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
        "xtick.color": MUTED, "ytick.color": MUTED,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "axes.grid.axis": "y", "grid.color": GRID, "grid.linewidth": 0.6,
        "legend.frameon": False,
    })


def save(fig, path: Path):
    fig.savefig(path)
    plt.close(fig)


def fraction_boxplots(comp: pd.DataFrame, cell_types: list[str], group_col: str,
                      order: list[str], stats_table: pd.DataFrame, title: str):
    """Per-sample cell-type fractions by group, one point per sample, q-value above each pair."""
    long = comp.melt(id_vars=[group_col], value_vars=cell_types, var_name="cell_type",
                     value_name="fraction")
    fig, ax = plt.subplots(figsize=(max(6, 0.75 * len(cell_types)), 4))
    common = dict(data=long, x="cell_type", y="fraction", hue=group_col, order=cell_types,
                  hue_order=order, palette=PAIR, ax=ax)
    sns.boxplot(**common, showfliers=False, width=0.7, linewidth=0.8, boxprops={"alpha": 0.35})
    sns.stripplot(**common, dodge=True, size=4, jitter=0.15, linewidth=0.5, edgecolor="white",
                  legend=False)
    top = long["fraction"].max()
    for i, ct in enumerate(cell_types):
        q = stats_table.loc[ct, "qvalue"]
        ax.text(i, top * 1.04, f"q={q:.2g}", ha="center", va="bottom", fontsize=7,
                color=INK if q < 0.1 else MUTED)
    ax.set_ylim(0, top * 1.12)
    ax.set_xlabel("")
    ax.set_ylabel("Fraction of sample's cells")
    ax.set_title(title, loc="left")
    ax.tick_params(axis="x", rotation=40)
    plt.setp(ax.get_xticklabels(), ha="right")
    ax.legend(title=None, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    return fig


def paired_slopes(long: pd.DataFrame, cell_types: list[str], stats_table: pd.DataFrame,
                  before: str = "Pre", after: str = "Post", time_col: str = "timepoint"):
    """Small multiples: each line is one patient's fraction before and after treatment."""
    n = len(cell_types)
    ncols = min(5, n)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(2.1 * ncols, 2.2 * nrows), squeeze=False)
    for ax, ct in zip(axes.flat, cell_types):
        wide = long.pivot(index="patient", columns=time_col, values=ct)
        for _, row in wide.iterrows():
            up = row[after] > row[before]
            ax.plot([0, 1], [row[before], row[after]], color=PAIR[1] if up else PAIR[0],
                    lw=1.2, alpha=0.8, marker="o", ms=3)
        q = stats_table.loc[ct, "qvalue"]
        ax.set_title(f"{ct}\nq={q:.2g}", fontsize=8, loc="left")
        ax.set_xticks([0, 1], [before, after])
        ax.set_xlim(-0.3, 1.3)
    for ax in axes.flat[n:]:
        ax.set_visible(False)
    handles = [plt.Line2D([], [], color=PAIR[1], label="increased"),
               plt.Line2D([], [], color=PAIR[0], label="decreased")]
    fig.legend(handles=handles, loc="upper right", ncols=2)
    fig.supylabel("Fraction of sample's cells", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return fig


def volcano(de: pd.DataFrame, title: str, group_a: str, group_b: str, n_labels: int = 8):
    """log2 fold change vs -log10 p; the top genes on each side are labeled."""
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    y = -np.log10(de["pvalue"])
    sig = de["qvalue"] < 0.1
    ax.scatter(de["log2fc"][~sig], y[~sig], s=4, color=GRID, linewidths=0, label="q ≥ 0.1")
    for sign, color, name in [(1, PAIR[1], f"higher in {group_b}"), (-1, PAIR[0], f"higher in {group_a}")]:
        pick = sig & (np.sign(de["log2fc"]) == sign)
        ax.scatter(de["log2fc"][pick], y[pick], s=8, color=color, linewidths=0,
                   label=f"{name} (q < 0.1, n={pick.sum()})")
        for gene in de[np.sign(de["log2fc"]) == sign].head(n_labels).index:
            ax.annotate(gene, (de.loc[gene, "log2fc"], y[gene]), fontsize=6.5, color=INK,
                        xytext=(3, 1), textcoords="offset points")
    ax.axvline(0, color=MUTED, lw=0.6)
    ax.grid(axis="x")
    ax.set_xlabel(f"log2 fold change ({group_b} vs {group_a})")
    ax.set_ylabel("-log10 p-value")
    ax.set_title(title, loc="left")
    ax.legend(loc="upper left", fontsize=7, markerscale=2)
    return fig


def score_by_group(values: pd.Series, labels: pd.Series, order: list[str], ylabel: str, title: str):
    """One point per sample, grouped, with the median marked."""
    df = pd.DataFrame({"value": values, "group": labels.loc[values.index]}).dropna()
    fig, ax = plt.subplots(figsize=(3.2, 3.8))
    sns.stripplot(df, x="group", y="value", order=order, hue="group", hue_order=order,
                  palette=PAIR, size=6, jitter=0.12, linewidth=0.5, edgecolor="white",
                  legend=False, ax=ax)
    for i, g in enumerate(order):
        med = df.loc[df["group"] == g, "value"].median()
        ax.plot([i - 0.25, i + 0.25], [med, med], color=INK, lw=1.5)
    counts = df["group"].value_counts()
    ax.set_xticks(range(len(order)), [f"{g}\n(n={counts.get(g, 0)})" for g in order])
    ax.set_xlabel("")
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left", fontsize=9)
    return fig
