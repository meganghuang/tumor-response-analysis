# Committed results

Output of one run of both pipelines at commit `832124a`, on GEO GSE120575
(md5 `8bb26ab1e694c1396de3751695fa90e8` for the expression archive). 14,505 cells and 35,491
genes after QC; 24 Leiden clusters and 9 labelled populations in the scanpy path, 22 Louvain
clusters and 10 populations in the Seurat path.

`results/` is otherwise gitignored, because the full output is regenerable and some of it is
large. **This is a partial set**, committed so the headline numbers are readable from the
repository without a two-hour re-run. Everything here was produced by the run described
above; nothing has been regenerated or edited by hand.

## What is here

```
figures/
  umap_cell_types.png                     scanpy: cell types and Leiden clusters
  violin_qc.png                           QC metrics, including the corrected and the
                                          uncorrected mitochondrial fraction side by side
  response_composition_per_patient.png    responder vs non-responder, one value per patient
  response_composition_all_samples.png    the same per biopsy (patients repeat)
  response_composition_baseline.png       the same, baseline biopsies only
  response_tcf7_cd8_per_patient.png       the paper's TCF7+ CD8 predictor, per patient
  response_tcf7_cd8_per_biopsy.png        the same, per biopsy
  response_signature_cd8_all.png          learned signature, CD8 T cells, cross-validated
  response_de_cd8.png                     volcano, CD8 T cells, per-biopsy means
  treatment_composition_paired.png        Pre vs Post in the 11 paired patients
  concordance_python_vs_seurat.png        the two paths compared
  seurat/                                 the same UMAP, TCF7 and signature plots from R

tables/
  response_signature_summary.csv          every AUC, with its unit and n
  response_composition_per_patient.csv    per-patient composition test
  treatment_composition_paired.csv        paired Pre/Post test
  composition_by_sample.csv               per-sample cell-type fractions
  response_de_cd8.csv                     full differential expression, CD8 T cells
  concordance_*.csv                       summary, composition, markers, AUCs, confusion
  seurat/                                 signature summary, per-patient composition, and
                                          the cell-type marker table
```

## What is absent

A workspace cleanup removed the files below before they were committed. Each is regenerable
by re-running the pipeline; none of them is referenced by the README's reported numbers
except the marker tables, which are discussed but not quoted.

- scanpy: `markers_by_cell_type.csv`, `markers_by_cluster.csv`, `umap_samples.png`, the two
  dot plots, `response_de_all_immune.{png,csv}`, the `all_immune` and `no_conflict`
  signature plots, the four `response_signature_scores_*.csv`,
  `tcf7_cd8_fraction_by_sample.csv`, `response_composition_{all_samples,baseline}.csv`, and
  the remaining `treatment_composition_*.{png,csv}`
- Seurat: everything except the three figures and three tables listed above

To regenerate the complete set:

```bash
python scripts/01_download_data.py
python scripts/02_qc_preprocess.py      # first run parses the text matrix
python scripts/03_cluster_annotate.py
python scripts/04_treatment_changes.py
python scripts/05_response_signatures.py
Rscript -e 'Sys.setenv(TUMOR_RESPONSE_ROOT=getwd()); source("R/02_qc_preprocess.R")'
Rscript -e 'Sys.setenv(TUMOR_RESPONSE_ROOT=getwd()); source("R/03_cluster_annotate.R")'
Rscript -e 'Sys.setenv(TUMOR_RESPONSE_ROOT=getwd()); source("R/04_treatment_changes.R")'
Rscript -e 'Sys.setenv(TUMOR_RESPONSE_ROOT=getwd()); source("R/05_response_signatures.R")'
python scripts/06_compare_paths.py
```

## One caveat on the Seurat marker table

`tables/seurat/markers_by_cell_type.csv` comes from `rank_sum_markers` in `R/00_common.R`,
not from `FindAllMarkers` — the same one-vs-rest Wilcoxon statistic, computed by ranking each
gene once and reusing those ranks across groups. Seurat's own implementation has no presto
available in this environment and falls back to a per-gene `wilcox.test` loop that does not
finish in hours at this cell count. The column names follow Seurat's convention
(`avg_log2FC`, `pct.1`, `pct.2`), with `z` and a Benjamini-Hochberg `p_val_adj` added.
