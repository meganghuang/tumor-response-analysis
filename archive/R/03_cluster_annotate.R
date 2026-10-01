## Clustering, marker-based cell-type labels, and data-driven marker genes.
##
## The marker sets and the labelling rule are the same as annotation.py: each cluster takes
## the marker set it scores highest on, unless the best score is below `min_score` or beats
## the runner-up by less than `min_margin`, in which case it is "Unresolved".
##
## Two implementation differences from the scanpy path, both reported in the concordance
## table rather than papered over: clusters come from Louvain (Seurat's default) where
## scanpy used Leiden, and module scores come from AddModuleScore where scanpy used
## score_genes. Both use control-gene bins, but the binning differs.

source(file.path("R", "00_common.R"))
suppressPackageStartupMessages(library(Seurat))

MARKER_GENES <- list(
  "B cells" = c("CD19", "MS4A1", "CD79A", "CD79B"),
  "Plasma cells" = c("IGJ", "JCHAIN", "MZB1", "SDC1"),
  "Monocytes-Macrophages" = c("CD14", "CD68", "LYZ", "C1QA", "C1QB", "CSF1R"),
  "Dendritic cells" = c("FCER1A", "CLEC10A", "CD1C"),
  "Plasmacytoid DCs" = c("LILRA4", "CLEC4C", "IL3RA"),
  "CD8 effector" = c("GZMH", "GNLY", "FGFBP2", "CX3CR1", "KLRG1"),
  "CD8 exhausted" = c("PDCD1", "HAVCR2", "LAG3", "ENTPD1", "TOX", "CXCL13"),
  "Memory-naive T" = c("TCF7", "IL7R", "CCR7", "SELL", "LEF1"),
  "Tregs" = c("FOXP3", "IL2RA", "CTLA4"),
  "NK cells" = c("NCAM1", "KLRF1", "SH2D1B"),
  "Cycling" = c("MKI67", "TOP2A", "TYMS")
)

NON_T_CELL_TYPES <- c("B cells", "Plasma cells", "Monocytes-Macrophages", "Dendritic cells",
                      "Plasmacytoid DCs", "NK cells")
UNRESOLVED <- "Unresolved"

#' Marker sets restricted to genes present in the object.
available_markers <- function(obj, markers = MARKER_GENES) {
  present <- lapply(markers, function(g) intersect(g, rownames(obj)))
  present[lengths(present) > 0]
}

#' Label each cluster by its top-scoring marker set.
score_cell_types <- function(obj, markers = MARKER_GENES, min_score = 0.2, min_margin = 0.1) {
  sets <- available_markers(obj, markers)
  for (ct in names(sets)) {
    obj <- AddModuleScore(obj, features = list(sets[[ct]]), name = "tmp_score", nbin = 24,
                          ctrl = 50, seed = 0)
    ## Carried across as a plain vector: Seurat matches metadata by name, and a vector
    ## named by anything other than cell barcodes is rejected as having no cell overlap.
    obj[[paste0("score_", ct)]] <- unname(obj[[]][["tmp_score1"]])
    obj[["tmp_score1"]] <- NULL
  }
  score_cols <- paste0("score_", names(sets))
  scores <- obj[[score_cols]]
  cluster_means <- aggregate(scores, by = list(cluster = obj$seurat_clusters), FUN = mean)
  mat <- as.matrix(cluster_means[, score_cols, drop = FALSE])
  best_idx <- apply(mat, 1, which.max)
  best <- mat[cbind(seq_len(nrow(mat)), best_idx)]
  second <- apply(mat, 1, function(r) sort(r, decreasing = TRUE)[2])
  labels <- names(sets)[best_idx]
  labels[best < min_score | (best - second) < min_margin] <- UNRESOLVED
  names(labels) <- as.character(cluster_means$cluster)
  obj$cell_type <- unname(labels[as.character(obj$seurat_clusters)])
  obj
}

#' CD8 T cells: CD8A or CD8B detected, outside the clusters that cannot be T cells.
cd8_mask <- function(obj, genes = c("CD8A", "CD8B")) {
  present <- intersect(genes, rownames(obj))
  if (!length(present)) stop("none of CD8A/CD8B are present in this annotation")
  expr <- LayerData(obj, assay = "RNA", layer = "data")[present, , drop = FALSE]
  detected <- Matrix::colSums(expr > 0) > 0
  setNames(detected & !(obj$cell_type %in% NON_T_CELL_TYPES), colnames(obj))
}

main <- function() {
  cfg <- load_cfg()
  processed_dir <- project_path(cfg$data$processed_dir)
  figures_dir <- r_figures_dir(cfg)
  tables_dir <- r_tables_dir(cfg)

  obj <- readRDS(file.path(processed_dir, "01_preprocessed_seurat.rds"))
  ## algorithm = 1 is Louvain; scanpy used Leiden. Recorded as a known difference.
  obj <- FindClusters(obj, resolution = cfg$clustering$resolution, algorithm = 1,
                      verbose = FALSE)
  obj <- score_cell_types(obj)
  print(table(obj$cell_type))

  png(file.path(figures_dir, "umap_cell_types.png"), width = 1800, height = 800, res = 150)
  print(DimPlot(obj, group.by = c("cell_type", "seurat_clusters"), label = TRUE,
                label.size = 2.5) & Seurat::NoAxes())
  dev.off()
  png(file.path(figures_dir, "umap_samples.png"), width = 1800, height = 700, res = 150)
  print(DimPlot(obj, group.by = c("response", "timepoint")) & Seurat::NoAxes())
  dev.off()
  png(file.path(figures_dir, "dotplot_known_markers.png"), width = 2200, height = 900, res = 150)
  print(DotPlot(obj, features = available_markers(obj), group.by = "seurat_clusters") +
          Seurat::RotatedAxis())
  dev.off()

  ## Seurat's conventional pre-filter (logfc.threshold = 0.25, min.pct = 0.1). Testing every
  ## gene, as the scanpy path does, is not tractable here: without presto the Wilcoxon test
  ## is pure R, and 35,000 genes x 24 clusters does not finish in hours. Genes that reach a
  ## top-50 marker list clear these thresholds comfortably, so the comparison is of the two
  ## top-50 sets rather than of the full rankings.
  ## One-vs-rest Wilcoxon over every gene, as the scanpy path does. `rank_sum_markers`
  ## computes the same statistic as FindAllMarkers(test.use = "wilcox") but ranks each gene
  ## once and reuses those ranks across groups; Seurat's own implementation falls back to a
  ## per-gene `wilcox.test` loop when presto is unavailable, which does not finish in hours
  ## at this cell count. Set TUMOR_RESPONSE_SEURAT_MARKERS=1 to use FindAllMarkers instead.
  m <- LayerData(obj, assay = "RNA", layer = "data")
  if (nzchar(Sys.getenv("TUMOR_RESPONSE_SEURAT_MARKERS"))) {
    Idents(obj) <- "seurat_clusters"
    by_cluster <- FindAllMarkers(obj, test.use = "wilcox", logfc.threshold = 0.25,
                                 min.pct = 0.1, verbose = FALSE)
    Idents(obj) <- "cell_type"
    by_type <- FindAllMarkers(obj, test.use = "wilcox", logfc.threshold = 0.25, min.pct = 0.1,
                              verbose = FALSE)
  } else {
    by_cluster <- rank_sum_markers(m, obj$seurat_clusters)
    by_type <- rank_sum_markers(m, obj$cell_type)
  }
  write.csv(do.call(rbind, by(by_cluster, by_cluster$cluster, head, 50)),
            file.path(tables_dir, "markers_by_cluster.csv"), row.names = FALSE)
  write.csv(do.call(rbind, by(by_type, by_type$cluster, head, 50)),
            file.path(tables_dir, "markers_by_cell_type.csv"), row.names = FALSE)

  ## Per-cell labels, so the two paths can be compared cell by cell.
  write.csv(data.frame(cell = colnames(obj), cluster = as.character(obj$seurat_clusters),
                       cell_type = obj$cell_type, sample = obj$sample),
            file.path(tables_dir, "cell_labels.csv"), row.names = FALSE)

  saveRDS(obj, file.path(processed_dir, "02_annotated_seurat.rds"), compress = FALSE)
}

main()
