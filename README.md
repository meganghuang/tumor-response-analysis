# Tumor Response Analysis

Practice replication of Sade-Feldman et al. 2018, *Cell*, "Defining T cell states associated
with response to checkpoint immunotherapy in melanoma" (GEO
[GSE120575](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120575)): 16,291 CD45+
immune cells from 48 melanoma biopsies of 32 patients treated with anti-PD1 and/or anti-CTLA4,
profiled with Smart-seq2.

Each sample is a baseline (`Pre`) or on-treatment (`Post`) biopsy labelled `Responder` or
`Non-responder`; 11 patients were biopsied at both timepoints.

## Questions

1. **Which cell populations increase or decrease with treatment?** (`04_treatment_changes.py`)
2. **What are the marker genes of each population?** (`03_cluster_annotate.py`)
3. **Which signatures distinguish responders from non-responders?** (`05_response_signatures.py`)

## Approach

- **Cell types.** Leiden clusters are labelled with the marker set they score highest on
  (`src/tumor_response/annotation.py`). Clusters with no clear winner are labelled
  `Unresolved` rather than forced into a weak match; in this dataset those are mixed CD4/CD8
  T-cell clusters.
- **One value per sample.** Cells from the same biopsy are not independent, so every group
  comparison uses per-sample values: cell-type fractions, or mean expression per sample
  ("pseudobulk"). Tests are Mann-Whitney U (unpaired) or Wilcoxon signed-rank (paired
  Pre/Post patients), with Benjamini-Hochberg q-values.
- **Honest signature evaluation.** Learned gene signatures are scored with
  leave-one-patient-out cross-validation: each patient's samples are scored by a signature
  derived without that patient, so the AUC is not inflated by fitting and testing on the
  same data.

## Results (current run)

**1. Treatment changes.** No population changes significantly between baseline and
on-treatment biopsies after multiple-testing correction, either across all samples or in the
11 paired patients. The strongest trend is a drop in NK cells in 9/11 paired patients
(p = 0.03, q = 0.35). This matches the paper's observation that composition tracks response
more than timepoint.

**2. Markers.** `results/tables/markers_by_cell_type.csv` and `markers_by_cluster.csv` list the
top 50 genes per population; `umap_cell_types.png` and the dot plots show them.

**3. Response signatures.**

| Signature | AUC |
|---|---|
| TCF7+ fraction of CD8 T cells (the paper's predictor, no fitting) | 0.72 |
| Learned 25+25 gene signature, CD8 T cells (cross-validated) | 0.85 |
| Learned 25+25 gene signature, all immune cells (cross-validated) | 0.89 |

- Responders have more memory/naive-like T cells (TCF7, IL7R, CCR7; q = 0.002) and B cells
  (q = 0.01).
- Non-responders have more cycling (q < 0.001) and exhausted CD8 T cells, monocytes/macrophages
  (q = 0.02), and plasmacytoid DCs (q = 0.04). The same pattern holds using baseline biopsies
  only.
- Within CD8 T cells, genes higher in responders include IL7R, CCR7, FOXP1 and ZBTB20
  (memory/stem-like); genes higher in non-responders include NKG7, HLA-DPA1 and the
  interferon/proteasome genes PSME1/PSME2 (activation/exhaustion).

**Caveats.** With 48 samples these are small-sample statistics. Some top DE hits are
pseudogenes or lncRNAs (e.g. `RP11-…`, `GAPDHP44`, `MALAT1`), which can reflect per-sample
technical effects as well as biology. Response is labelled per biopsy, and a few patients have
several biopsies, so samples are not fully independent (the cross-validation holds out whole
patients to account for this).

## Project layout

```
config/config.yaml        # dataset URLs, QC thresholds, analysis parameters
data/raw/                 # downloaded GEO files (not tracked by git)
data/processed/           # intermediate .h5ad files (not tracked by git)
notebooks/                # exploratory analysis
scripts/                  # numbered pipeline steps, run in order
src/tumor_response/       # reusable analysis code
  io.py                   #   GSE120575 loader
  qc.py, preprocess.py    #   QC filtering, PCA / neighbors / UMAP
  annotation.py           #   clustering, cell-type labels, marker genes, CD8 selection
  composition.py          #   per-sample cell-type fractions and group tests
  signatures.py           #   pseudobulk DE, signature scoring, cross-validation
  plotting.py             #   figures
results/figures/          # generated plots (not tracked by git)
results/tables/           # generated tables (not tracked by git)
tests/                    # unit tests on synthetic data
```

## Setup

```bash
python -m venv .venv
.venv/Scripts/activate      # Windows; use `source .venv/bin/activate` on macOS/Linux
pip install -e ".[dev]"
```

## Running the pipeline

```bash
python scripts/01_download_data.py      # ~120 MB from GEO
python scripts/02_qc_preprocess.py      # first run parses the text matrix (~10 min), then cached
python scripts/03_cluster_annotate.py
python scripts/04_treatment_changes.py
python scripts/05_response_signatures.py
```

Each step reads the previous step's output from `data/processed/` and writes figures and
tables to `results/`.

## Tests

```bash
pytest
```
