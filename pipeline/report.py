
"""Assemble results/summary_report.md from the computed result objects.

Every number in the report is formatted from a DataFrame produced by the
pipeline; nothing is hard-coded, so the report cannot drift from the tables.
"""
from __future__ import annotations

import datetime as _dt
import pathlib

import pandas as pd


def md_table(df: pd.DataFrame, floatfmt: str = "{:.3g}", max_rows: int | None = None) -> str:
    if max_rows is not None:
        df = df.head(max_rows)
    def fmt(v):
        if isinstance(v, float):
            return floatfmt.format(v)
        return str(v)
    head = "| " + " | ".join(str(c) for c in df.columns) + " |"
    rule = "| " + " | ".join("---" for _ in df.columns) + " |"
    rows = ["| " + " | ".join(fmt(v) for v in r) + " |" for r in df.itertuples(index=False)]
    return "\n".join([head, rule] + rows)


def build(res: dict, out: pathlib.Path) -> pathlib.Path:
    d = res
    L: list[str] = []
    A = L.append

    A("# Cell-population dynamics and response signatures in melanoma checkpoint immunotherapy")
    A("")
    A(f"**Dataset** GSE120575 — Sade-Feldman *et al.*, *Cell* 175:998-1013 (2018). "
      f"CD45+ cells from melanoma biopsies, Smart-seq2, distributed as TPM.")
    A(f"**Pipeline** `pipeline/` v{d['version']} · generated "
      f"{_dt.datetime.now(_dt.timezone.utc):%Y-%m-%d} · random seed "
      f"{d['cfg']['embedding']['random_state']}")
    A("")
    A("## Headline findings")
    A("")
    for b in d["headline"]:
        A(f"- {b}")
    A("")

    A("## 1. Data and quality control")
    A("")
    A(f"The GEO matrix contains {d['n_cells_raw']:,} cells x {d['n_genes_raw']:,} genes from "
      f"{d['n_samples_raw']} biopsies of {d['n_patients_raw']} patients. "
      f"After filtering, {d['n_cells']:,} cells ({100*d['n_cells']/d['n_cells_raw']:.1f}%) "
      f"and {d['n_genes']:,} genes from {d['n_samples']} biopsies of {d['n_patients']} patients "
      f"enter the analysis.")
    A("")
    A("Thresholds (`config/config.yaml`): "
      + ", ".join(f"`{k}` = {v}" for k, v in d["cfg"]["qc"].items()) + ".")
    A("")
    A("![QC distributions and the cell count surviving each filter](qc_metrics.png)")
    A("")
    A(md_table(d["qc_summary"].drop(columns=["n_genes_retained"])))
    A("")

    A("## 2. Clustering and annotation")
    A("")
    A(f"Log2(TPM/10+1) expression, {d['cfg']['embedding']['n_top_genes']} highly variable genes, "
      f"{d['cfg']['embedding']['n_pcs']} principal components, a "
      f"{d['cfg']['embedding']['n_neighbors']}-nearest-neighbour graph and Leiden clustering at "
      f"resolution {d['cfg']['embedding']['leiden_resolution']} give "
      f"{d['n_clusters']} clusters.")
    A("")
    A("![UMAP coloured by cluster, timepoint and response](umap_clusters.png)")
    A("")
    A(f"Clusters are labelled by the argmax of z-scored canonical-programme scores. "
      f"{d['pct_unresolved']:.1f}% of cells fall in clusters where no programme scores above "
      f"the cross-cluster mean and are left Unresolved. CD8 T cells "
      f"({d['n_cd8']:,} cells) were re-embedded and sub-clustered at resolution "
      f"{d['cfg']['embedding']['cd8_leiden_resolution']} to resolve differentiation states.")
    A("")
    A("![Annotated lineages and CD8 T-cell states](umap_celltypes.png)")
    A("")
    A(md_table(d["population_counts"]))
    A("")

    A("## 3. Marker genes")
    A("")
    A("![Canonical lineage markers per population](marker_dotplot.png)")
    A("")
    A("Data-driven markers come from a one-vs-rest Wilcoxon rank-sum test per population. "
      f"All {d['n_marker_rows']:,} gene-population pairs were tested; `cluster_markers.csv` "
      f"holds the {d['n_marker_sig']:,} that are enriched in their population at BH q < 0.05 "
      f"and detected in at least 10% of its cells, covering {d['n_marker_pops']} "
      f"of {len(d['population_counts'])} populations.")
    A("")
    A("![Top data-driven markers, z-scored mean expression](marker_heatmap.png)")
    A("")

    A("## 4. Populations that expand or contract with treatment")
    A("")
    A(f"Paired Wilcoxon signed-rank on population fractions, pre- versus on-treatment, "
      f"restricted to the {d['n_paired_patients']} patients who contributed both. "
      f"The paired design is the only way to read a treatment effect here: biopsy timing is "
      f"confounded with patient.")
    A("")
    A("![Composition by response and timepoint, and per biopsy](composition_barplots.png)")
    A("")
    A(md_table(d["treatment_table"]))
    A("")
    A(d["treatment_prose"])
    A("")

    A("## 5. Populations associated with response")
    A("")
    A("The patient is the unit of analysis. Biopsies from one patient are not independent, "
      "so each patient contributes one fraction per population (the mean over their biopsies) "
      "and groups are compared with Mann-Whitney U, BH-corrected across populations.")
    A("")
    A("![Per-population abundance, responders vs non-responders](composition_boxplots.png)")
    A("")
    A(md_table(d["response_table"]))
    A("")
    A(d["response_prose"])
    A("")

    A("## 6. Differential expression, responder vs non-responder")
    A("")
    A(f"Patient-level pseudobulk (mean log2(TPM/10+1) over a patient's cells), Welch t-test, "
      f"BH across genes. Full tables: `responder_DE_all.csv` "
      f"({d['n_de_all']:,} genes, {d['n_de_all_sig']:,} at q < "
      f"{d['cfg']['differential_expression']['fdr_alpha']}) and `responder_DE_CD8.csv` "
      f"({d['n_de_cd8']:,} genes, {d['n_de_cd8_sig']:,} at q < "
      f"{d['cfg']['differential_expression']['fdr_alpha']}).")
    A("")
    A("**Top genes, all immune cells**")
    A("")
    A(md_table(d["de_all_top"]))
    A("")
    A("**Top genes, CD8 T cells**")
    A("")
    A(md_table(d["de_cd8_top"]))
    A("")

    A("## 7. Signatures that stratify responders")
    A("")
    A("A signature fitted on all patients and scored on those same patients is circular. "
      "Every AUC below is leave-one-patient-out: gene selection and gene-wise standardisation "
      "are re-derived from the training patients alone and the held-out patient is scored "
      "with them. P-values come from "
      f"{d['cfg']['signatures']['n_permutations']} label permutations.")
    A("")
    A("![ROC curves, leave-one-patient-out](signature_auc.png)")
    A("")
    A(md_table(d["auc_table"]))
    A("")
    A("![Signature score distributions](signature_scores.png)")
    A("")
    A(d["signature_prose"])
    A("")
    A("Final gene lists (fitted on all patients, for scoring a future cohort) are in "
      "`signature_genes.json`.")
    A("")

    A("## 8. Methods summary")
    A("")
    for m in d["methods"]:
        A(f"- {m}")
    A("")

    A("## 9. Limitations")
    A("")
    for m in d["limitations"]:
        A(f"- {m}")
    A("")

    A("## 10. Reproducing this report")
    A("")
    A("```bash")
    A("conda env create -f environment.yml && conda activate melanoma-icb")
    A("python run_pipeline.py            # downloads GSE120575 and writes results/")
    A("```")
    A("")
    A("Artifacts written to `results/`: "
      + ", ".join(f"`{f}`" for f in d["artifacts"]) + ".")
    A("")

    out.write_text("\n".join(L), encoding="utf-8")
    return out
