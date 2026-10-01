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

## Two implementations

The analysis is implemented twice against the same raw files, the same thresholds and the
same statistics: in Python/scanpy under `src/tumor_response/` and `scripts/`, and in R/Seurat
under `R/`. `scripts/06_compare_paths.py` measures how far the two land apart. See
[`R/README.md`](R/README.md) for the Seurat path and the places the two legitimately differ.

## Corrections

`AUDIT.md` lists every defect found in the first version of this pipeline, what it did to the
reported numbers, and how it was fixed. The four that changed results:

1. **The mitochondrial fraction was computed from log-scale values**, making it a ratio of
   sums of `ln(TPM+1)` rather than a fraction of transcripts. Corrected, the median cell is
   12.7% mitochondrial rather than an apparent 1%, and `max_pct_mito: 20` removes 1,786 cells
   (11%) where before it removed none. Both implementations agree on this to the cell.
2. **Forty-eight biopsies from thirty-two patients were treated as independent.** Every
   response comparison is now reported per patient (primary) as well as per biopsy.
3. **Four patients carry contradictory response labels** across their lesions (P1, P4, P5,
   P28) and so appeared in both arms of every comparison. They are now excluded from the
   patient-level analysis as `Conflicting`, and the cross-validation is reported with and
   without them.
4. **The gene detection filter sat outside the cross-validation fold.** It is now refitted on
   training samples only, alongside the gene selection and the z-reference.

Neither gene-count filter does anything on this dataset: detected genes run from 1,093 to
9,315, so `min_genes: 1000` and `max_genes: 10000` remove no cells.

## scanpy versus Seurat

Both paths keep exactly the same 14,505 cells: the two loaders agree on 16,291 cells, 37
mitochondrial genes, a mean mitochondrial fraction of 13.94%, and the same 1,786 cells
removed. Downstream, they agree on the answers and disagree on the labels.

| Comparison | Value |
|---|---|
| Cells in both paths | 14,505 |
| Clusters | 24 (Leiden) vs 22 (Louvain) |
| Adjusted Rand index, clusters | 0.29 |
| Adjusted mutual information, clusters | 0.50 |
| Cells given the same cell-type label | 47.5% |
| Median per-sample composition Pearson r | 0.985 |
| Median top-50 marker Jaccard | 0.27 |
| Largest response-AUC difference | 0.008 |

Every response signature lands within 0.008 AUC of the other path, including the paper's
TCF7+ CD8 predictor (0.789 vs 0.789 per patient). Cell-type *abundances* agree closely for
populations defined by unambiguous lineage markers — monocytes/macrophages r = 0.999, pDCs
0.997, plasma cells 0.995, B cells and cycling 0.985 and 0.988 — and poorly for the T-cell
states: CD8 effector r = 0.64, Tregs 0.60, memory/naive 0.76.

The reason is worth acting on. **The Seurat path assigns 3,322 cells (23%) to a population the
scanpy path does not label at all**: 44% of the cells scanpy calls `Unresolved` are labelled
`CD8 exhausted` by Seurat, and `Unresolved` falls from 45.9% of cells to 4.9%. Under Leiden
clustering with `min_margin: 0.1`, no scanpy cluster wins that marker set cleanly. The
labelling rule, not the data, is the limiting step.

The label itself should be read with care rather than taken at face value. That group's own
top differential markers are `CD8A`, `CCL5`, `CD8B`, `NKG7` and `KLRK1` — generic CD8 and
cytotoxic genes, not exhaustion-specific ones. What the module score identifies is a large
cytotoxic CD8 compartment that scores highest on the exhaustion set, not a population
independently confirmed as exhausted. Distinguishing exhausted from effector CD8 states here
needs a better-separated clustering or a graded score, not a winner-takes-all label.

That the response AUCs survive this is informative: they are driven by pseudobulk expression
and by a marker-based CD8 definition, neither of which depends on the cluster labels.

## Results (current run)

14,505 cells and 35,491 genes survive QC, in 24 Leiden clusters and 9 labelled populations.

**1. Treatment changes.** No population changes significantly between baseline and
on-treatment biopsies. In the 11 paired patients every population returns q = 1.0; restricted
to the 21 patients biopsied at a single timepoint, the strongest signals are an increase in
monocytes/macrophages and CD8 effector cells (both q = 0.13). This matches the paper's
observation that composition tracks response more than timepoint. Two patients, P1 and P6,
were biopsied under more than one therapy, so their Pre/Post pair is not one treatment course.

**2. Markers.** `results/tables/markers_by_cell_type.csv` and `markers_by_cluster.csv` list the
top 50 genes per population; `umap_cell_types.png` and the dot plots show them. These p-values
come from a cell-level test and are descriptive, not inferential — see AUDIT.md section B3.

**3. Response signatures.**

| Signature | Unit | n | AUC |
|---|---|---|---|
| TCF7+ fraction of CD8 T cells (the paper's predictor, no fitting) | patient | 28 | 0.79 |
| TCF7+ fraction of CD8 T cells | biopsy | 47 | 0.72 |
| Learned 25+25 gene signature, all immune cells (cross-validated) | biopsy | 48 | 0.88 |
| Learned 25+25 gene signature, all immune cells, conflicted patients dropped | biopsy | 38 | 0.89 |
| Learned 25+25 gene signature, CD8 T cells (cross-validated) | biopsy | 47 | 0.89 |
| Learned 25+25 gene signature, CD8 T cells, conflicted patients dropped | biopsy | 38 | 0.89 |

Moving to the patient as the unit *strengthens* the paper's predictor (0.72 to 0.79) and the
composition differences, despite cutting the number of observations from 48 to 28.

Per patient, five populations differ at q < 0.1:

| Population | Higher in | q |
|---|---|---|
| Cycling | non-responders | 0.0015 |
| B cells | responders | 0.0076 |
| Memory-naive T | responders | 0.018 |
| CD8 effector | non-responders | 0.027 |
| Plasmacytoid DCs | non-responders | 0.043 |

Within CD8 T cells, genes higher in responders include IL7R, CCR7, FOXP1 and ZBTB20
(memory/stem-like); genes higher in non-responders include NKG7, HLA-DPA1 and the
interferon/proteasome genes PSME1/PSME2 (activation/exhaustion).

**Caveats.** With 32 patients these are small-sample statistics. 45% of cells (6,572) fall in
`Unresolved`, so the labelled fractions describe a little over half the data; no cluster wins
the NK or exhausted-CD8 marker sets under the current margin rule, so those populations are
absent from the composition tables entirely. Some top differential-expression hits are
pseudogenes or lncRNAs (`RP11-…`, `GAPDHP44`, `MALAT1`), which can reflect per-sample technical
effects as well as biology. Cell-type fractions are compositional and are tested one at a
time, which ignores the constraint that they sum to one.

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
