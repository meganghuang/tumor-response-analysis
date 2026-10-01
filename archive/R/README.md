# Seurat path

A parallel implementation of the same analysis in R/Seurat, reading the same raw GEO files
and applying the same QC thresholds, marker sets, statistics and cross-validation scheme as
`src/tumor_response/`. Its outputs go to `results/figures/seurat/` and
`results/tables/seurat/` so neither path overwrites the other, and
`scripts/06_compare_paths.py` measures how far apart they land.

## Files

```
R/00_common.R              config, paths, and every statistic shared by the steps below
R/01_load.R                GSE120575 loader: decompress once, read gene blocks, sparsify
R/02_qc_preprocess.R       QC filtering, Seurat object, variable genes, PCA, neighbours, UMAP
R/03_cluster_annotate.R    clustering, marker-based labels, marker genes, per-cell label export
R/04_treatment_changes.R   goal 1: composition change between baseline and on treatment
R/05_response_signatures.R goal 3: responder signatures, TCF7+ CD8, cross-validated scores
```

## Running

```r
Sys.setenv(TUMOR_RESPONSE_ROOT = getwd())   # run from the project root
source("R/02_qc_preprocess.R")              # sources 01_load.R; first run parses the matrix
source("R/03_cluster_annotate.R")
source("R/04_treatment_changes.R")
source("R/05_response_signatures.R")
```

Step 02 decompresses the expression archive to `data/raw/*.txt` on first run (about 1.4 GB,
gitignored) and caches the parsed matrix as `data/processed/00_raw_seurat.rds`.

## Loading notes

The expression matrix is ~55,700 genes x 16,291 cells, which needs roughly 7 GB dense, so
`01_load.R` reads 4,000 gene rows at a time and converts each block to a sparse matrix before
moving on. Values arrive as log2(TPM+1) and are multiplied by ln(2) to give ln(TPM+1), the
scale Seurat's `data` layer expects. They are **not** counts: `NormalizeData` must not be run,
and the values are written into the `data` layer directly with `LayerData()<-`.

## Where the two paths legitimately differ

These are implementation differences, not disagreements about method. They are measured in
`results/tables/concordance_summary.csv` rather than assumed away.

| Step | scanpy | Seurat |
|---|---|---|
| Variable genes | `highly_variable_genes` flavor `seurat` (dispersion bins) | `FindVariableFeatures` `vst` |
| Community detection | Leiden (`flavor="igraph"`, 2 iterations) | Louvain (`algorithm = 1`) |
| Neighbour graph | exact kNN via scikit-learn | `FindNeighbors` (Annoy, approximate) |
| UMAP | `scanpy.tl.umap` | `RunUMAP` (uwot) |
| Module scores | `score_genes` | `AddModuleScore` |
| Marker genes | Wilcoxon on all genes | Wilcoxon with `min.pct = 0.1` |

Leiden is not used on the R side because `FindClusters(algorithm = 4)` needs `leidenAlg`,
which is not installed here; the cluster-level agreement between Leiden and Louvain is
reported as an adjusted Rand index instead of being hidden.

## Environment

Seurat 5.5.1 under R 4.4.3, in a conda environment created with `r-base=4.4`, `r-seurat`,
`r-matrix`, `r-yaml` and `r-data.table` from conda-forge.

Two things that do not work on this platform, recorded so they are not retried:

- **conda-forge `r-seurat` alongside `r-tidyverse`** produces an R that dies at startup with
  a Mingw-w64 pseudo-relocation failure. Installing Seurat without tidyverse in the same
  environment is fine.
- **CRAN binaries installed into a user library** install cleanly but cannot be loaded: the
  sandbox refuses to `dyn.load` a shared library from a writable path, so every compiled
  dependency fails. Seurat has to come from conda, not `install.packages()`.

`R.utils` is absent, so `fread()` cannot read `.gz` directly; `01_load.R` decompresses with
base R instead, which is also faster when reading the file in blocks.
