# Melanoma checkpoint-immunotherapy scRNA-seq pipeline

A reproducible, from-scratch analysis of **GSE120575** (Sade-Feldman *et al.*,
*Cell* 175:998-1013, 2018): 16,291 CD45+ cells profiled by Smart-seq2 from 48
melanoma biopsies of 32 patients treated with checkpoint blockade.

The pipeline answers three questions:

1. **Which immune populations expand or contract with treatment?**
   Paired pre- versus on-treatment comparison within each patient.
2. **What defines each population?**
   One-vs-rest marker genes alongside canonical lineage panels.
3. **Which transcriptional signatures stratify responders?**
   Patient-level pseudobulk differential expression, with every reported AUC
   estimated leave-one-patient-out.

## Quick start

```bash
conda env create -f environment.yml
conda activate melanoma-icb
python run_pipeline.py
```

The first run downloads ~127 MB from GEO and parses the TPM text matrix
(~10 min); the parsed matrix is cached to `data/processed/raw_tpm.h5ad`, so
later runs skip straight to the analysis. Everything lands in `results/`.

## Layout

```
config/config.yaml     every threshold and parameter, in one place
pipeline/              the library
  data.py              GEO download, matrix/metadata parsing
  qc.py                per-cell quality control
  embed.py             normalisation, HVGs, PCA, kNN, Leiden, UMAP
  annotate.py          lineage programmes and CD8 differentiation states
  markers.py           one-vs-rest marker genes
  composition.py       treatment (paired) and response (patient-level) tests
  de.py                patient pseudobulk differential expression
  signatures.py        signature learning and leave-one-patient-out evaluation
  figures.py           all figures
  report.py            assembles summary_report.md from the result objects
run_pipeline.py        end-to-end driver
tests/                 unit tests on synthetic data
results/               generated artifacts (committed)
archive/               the previous iteration of this project, kept for reference
```

## Outputs

| File | Contents |
| --- | --- |
| `results/summary_report.md` | Full write-up with every figure embedded |
| `results/qc_metrics.png` | QC distributions and the filter cascade |
| `results/umap_clusters.png` | UMAP by cluster, timepoint and response |
| `results/umap_celltypes.png` | Annotated lineages and CD8 states |
| `results/marker_dotplot.png` | Canonical markers per population |
| `results/marker_heatmap.png` | Top data-driven markers |
| `results/composition_barplots.png` | Stacked composition by timepoint and response |
| `results/composition_boxplots.png` | Per-population abundance, R vs NR, with p-values |
| `results/signature_auc.png` | Leave-one-patient-out ROC curves |
| `results/signature_scores.png` | Signature score distributions |
| `results/composition_stats.csv` | Treatment and response composition statistics |
| `results/cluster_markers.csv` | Per-population marker genes |
| `results/responder_DE_all.csv` | Pseudobulk DE, all immune cells |
| `results/responder_DE_CD8.csv` | Pseudobulk DE, CD8 T cells |
| `results/signature_genes.json` | Final gene lists and AUCs |

## Statistical stance

Three choices drive most of the difference between this pipeline and a naive one:

- **The patient is the unit of analysis** for every response comparison. Cells
  within a biopsy, and biopsies within a patient, are not independent
  replicates; testing at the cell level inflates significance by two orders of
  magnitude without adding information.
- **Treatment effects are paired.** Timepoint is confounded with patient, so
  pre/on-treatment comparisons are made within patient and tested with a
  signed-rank test.
- **Signature performance is leave-one-patient-out.** Gene selection and
  standardisation are re-derived from the training patients at every fold, and
  p-values come from label permutation.

## Tests

```bash
pytest -q
```

The tests run on synthetic data and check the statistics and annotation logic
(filter cascade, log transform, Unresolved rule, paired and unpaired tests,
pseudobulk collapse, and that leave-one-patient-out never sees the held-out
patient's label).
