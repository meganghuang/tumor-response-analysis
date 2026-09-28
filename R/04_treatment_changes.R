## Goal 1: cell populations that increase or decrease with treatment (Pre vs Post).
##
## The primary unpaired test uses only patients biopsied at a single timepoint, so the two
## groups share no patient; the paired test uses the patients biopsied at both. The
## all-samples comparison is kept for continuity with the original run and reported as
## non-independent.

source(file.path("R", "00_common.R"))
suppressPackageStartupMessages(library(Seurat))

main <- function() {
  cfg <- load_cfg()
  processed_dir <- project_path(cfg$data$processed_dir)
  figures_dir <- r_figures_dir(cfg)
  tables_dir <- r_tables_dir(cfg)

  obj <- readRDS(file.path(processed_dir, "02_annotated_seurat.rds"))
  meta <- obj[[]]
  comp <- composition_by_sample(meta)
  write.csv(comp, file.path(tables_dir, "composition_by_sample.csv"), row.names = FALSE)
  cell_types <- sort(unique(meta$cell_type))

  disjoint <- single_timepoint_samples(comp)
  independent <- compare_fractions(disjoint, cell_types, "timepoint", "Pre", "Post")
  write.csv(independent, file.path(tables_dir, "treatment_composition_independent_patients.csv"),
            row.names = FALSE)
  n <- attr(independent, "n")
  cat(sprintf("Independent patients only: Pre n=%d, Post n=%d, %d cell types at q < 0.1\n",
              n[["Pre"]], n[["Post"]], sum(independent$qvalue < 0.1, na.rm = TRUE)))

  unpaired <- compare_fractions(comp, cell_types, "timepoint", "Pre", "Post")
  write.csv(unpaired, file.path(tables_dir, "treatment_composition_all_samples.csv"),
            row.names = FALSE)

  paired <- paired_fractions(comp, cell_types)
  write.csv(paired, file.path(tables_dir, "treatment_composition_paired.csv"), row.names = FALSE)
  cat(sprintf("Paired patients: %d\n", attr(paired, "n_patients")))
  switchers <- inconsistent_patients(comp, "therapy")
  if (length(switchers)) {
    cat("Patients whose biopsies span more than one therapy, so their Pre/Post pair is not",
        "one treatment course:", paste(switchers, collapse = ", "), "\n")
  }
  print(paired, digits = 3)

  ## Paired slopes: one line per patient per cell type.
  pre <- attr(paired, "pre"); post <- attr(paired, "post")
  png(file.path(figures_dir, "treatment_composition_paired.png"),
      width = 400 * min(5, length(cell_types)),
      height = 420 * ceiling(length(cell_types) / 5), res = 150)
  op <- par(mfrow = c(ceiling(length(cell_types) / 5), min(5, length(cell_types))),
            mar = c(3, 4, 3, 1))
  for (ct in cell_types) {
    rng <- range(c(pre[[ct]], post[[ct]]))
    plot(NA, xlim = c(0.8, 2.2), ylim = rng, xaxt = "n", xlab = "",
         ylab = "Fraction of sample's cells",
         main = sprintf("%s\nq=%.2g", ct, paired$qvalue[paired$cell_type == ct]),
         cex.main = 0.9)
    axis(1, at = c(1, 2), labels = c("Pre", "Post"))
    up <- post[[ct]] > pre[[ct]]
    segments(1, pre[[ct]], 2, post[[ct]], col = ifelse(up, "#eb6834", "#2a78d6"), lwd = 1.4)
    points(rep(1, nrow(pre)), pre[[ct]], pch = 16, cex = 0.6, col = "#6b6a64")
    points(rep(2, nrow(post)), post[[ct]], pch = 16, cex = 0.6, col = "#6b6a64")
  }
  par(op)
  dev.off()
}

main()
