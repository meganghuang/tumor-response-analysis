## QC filtering and dimensionality reduction with Seurat, at the parameters in config.yaml.
##
## The mitochondrial fraction is computed on the linear (TPM) scale, matching the corrected
## Python path: summing ln(TPM+1) would give a ratio of logarithms, not a fraction of
## transcripts. Values are placed in the `data` layer directly -- they are TPM, not counts,
## so NormalizeData must not be run.

source(file.path("R", "00_common.R"))
source(file.path("R", "01_load.R"))
suppressPackageStartupMessages(library(Seurat))

#' Percent of TPM contributed by mitochondrial genes, computed on the linear scale.
mito_percent <- function(m) {
  mt <- grepl("^MT-", toupper(rownames(m)))
  linear <- m
  linear@x <- expm1(linear@x)
  total <- Matrix::colSums(linear)
  mito <- if (any(mt)) Matrix::colSums(linear[mt, , drop = FALSE]) else rep(0, ncol(linear))
  list(pct = ifelse(total > 0, 100 * mito / total, 0), n_mito_genes = sum(mt))
}

#' The uncorrected log-scale ratio, kept only so the two can be compared.
mito_percent_logscale <- function(m) {
  mt <- grepl("^MT-", toupper(rownames(m)))
  total <- Matrix::colSums(m)
  mito <- if (any(mt)) Matrix::colSums(m[mt, , drop = FALSE]) else rep(0, ncol(m))
  ifelse(total > 0, 100 * mito / total, 0)
}

main <- function() {
  cfg <- load_cfg()
  processed_dir <- project_path(cfg$data$processed_dir)
  figures_dir <- r_figures_dir(cfg)
  raw_dir <- project_path(cfg$data$raw_dir)

  cache <- file.path(processed_dir, "00_raw_seurat.rds")
  if (file.exists(cache)) {
    data <- readRDS(cache)
  } else {
    data <- load_sade_feldman(raw_dir)
    saveRDS(data, cache, compress = FALSE)
  }
  m <- data$matrix
  meta <- data$meta
  cat(sprintf("Loaded: %d cells x %d genes, %d samples\n", ncol(m), nrow(m),
              length(unique(meta$sample))))

  mito <- mito_percent(m)
  meta$pct_counts_mt <- mito$pct
  meta$pct_counts_mt_logscale <- mito_percent_logscale(m)
  meta$n_genes_by_counts <- Matrix::colSums(m > 0)
  cat(sprintf("Mitochondrial genes in this annotation: %d\n", mito$n_mito_genes))
  cat(sprintf("pct mito, linear scale vs uncorrected log scale (mean, 95th pct): %.2f / %.2f vs %.2f / %.2f\n",
              mean(meta$pct_counts_mt), quantile(meta$pct_counts_mt, 0.95),
              mean(meta$pct_counts_mt_logscale), quantile(meta$pct_counts_mt_logscale, 0.95)))

  q <- cfg$qc
  removed <- c(too_few_genes = sum(meta$n_genes_by_counts < q$min_genes),
               too_many_genes = sum(meta$n_genes_by_counts > q$max_genes),
               high_mito = sum(meta$pct_counts_mt > q$max_pct_mito))
  keep <- meta$n_genes_by_counts >= q$min_genes & meta$n_genes_by_counts <= q$max_genes &
    meta$pct_counts_mt <= q$max_pct_mito
  cat("Cells removed by each QC rule:", paste(names(removed), removed, sep = "=", collapse = ", "),
      sprintf("| total=%d\n", sum(!keep)))
  m <- m[, keep, drop = FALSE]
  meta <- meta[keep, , drop = FALSE]
  gene_cells <- Matrix::rowSums(m > 0)
  m <- m[gene_cells >= q$min_cells, , drop = FALSE]
  cat(sprintf("After QC: %d cells x %d genes\n", ncol(m), nrow(m)))

  ## Seurat rewrites underscores in feature names to dashes, which would leave the matrix
  ## assigned into the `data` layer below with names that no longer match the object.
  renamed <- make.unique(gsub("_", "-", rownames(m)))
  cat(sprintf("Gene names adjusted to Seurat's convention: %d\n", sum(renamed != rownames(m))))
  rownames(m) <- renamed

  obj <- CreateSeuratObject(counts = m, meta.data = meta, min.cells = 0, min.features = 0)
  ## The values are TPM already on a log scale, so they go straight into `data`.
  LayerData(obj, assay = "RNA", layer = "data") <- m
  p <- cfg$preprocess
  ## "dispersion" ranks genes on the `data` layer, which is the correct scale here; "vst"
  ## would model the `counts` layer, and these values are ln(TPM+1), not counts. This is the
  ## closest available match to scanpy's dispersion-based `seurat` flavor.
  obj <- FindVariableFeatures(obj, selection.method = "dispersion", nfeatures = p$n_top_genes,
                              verbose = FALSE)
  obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
  obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = p$n_pcs, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = seq_len(p$n_pcs), k.param = cfg$clustering$n_neighbors,
                       verbose = FALSE)
  obj <- RunUMAP(obj, dims = seq_len(p$n_pcs), verbose = FALSE)

  png(file.path(figures_dir, "violin_qc.png"), width = 1500, height = 550, res = 150)
  print(VlnPlot(obj, features = c("n_genes_by_counts", "pct_counts_mt",
                                  "pct_counts_mt_logscale"), pt.size = 0, ncol = 3) &
          Seurat::NoLegend())
  dev.off()

  saveRDS(obj, file.path(processed_dir, "01_preprocessed_seurat.rds"), compress = FALSE)
}

main()
