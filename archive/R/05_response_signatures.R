## Goal 3: signatures that distinguish responders from non-responders.
##
## Mirrors the corrected Python path: composition per patient (primary) and per biopsy
## (secondary), the TCF7+ CD8 fraction at both units, pseudobulk differential expression,
## and leave-one-patient-out cross-validation with the detection filter, gene selection and
## z-reference all refitted inside each fold.

source(file.path("R", "00_common.R"))
suppressPackageStartupMessages(library(Seurat))

NR <- "Non-responder"
R_ <- "Responder"

NON_T_CELL_TYPES <- c("B cells", "Plasma cells", "Monocytes-Macrophages", "Dendritic cells",
                      "Plasmacytoid DCs", "NK cells")

cd8_mask <- function(obj, genes = c("CD8A", "CD8B")) {
  present <- intersect(genes, rownames(obj))
  if (!length(present)) stop("none of CD8A/CD8B are present in this annotation")
  expr <- LayerData(obj, assay = "RNA", layer = "data")[present, , drop = FALSE]
  detected <- Matrix::colSums(expr > 0) > 0
  setNames(detected & !(obj$cell_type %in% NON_T_CELL_TYPES), colnames(obj))
}

main <- function() {
  cfg <- load_cfg()
  a <- cfg$analysis
  processed_dir <- project_path(cfg$data$processed_dir)
  figures_dir <- r_figures_dir(cfg)
  tables_dir <- r_tables_dir(cfg)

  obj <- readRDS(file.path(processed_dir, "02_annotated_seurat.rds"))
  meta <- obj[[]]
  m <- LayerData(obj, assay = "RNA", layer = "data")

  comp <- composition_by_sample(meta)
  cell_types <- sort(unique(meta$cell_type))
  sample_info <- unique(meta[, c("sample", "response", "patient", "timepoint")])
  rownames(sample_info) <- sample_info$sample
  labels <- setNames(sample_info$response, sample_info$sample)
  patients <- setNames(sample_info$patient, sample_info$sample)
  patient_labels <- vapply(sort(unique(sample_info$patient)), function(p)
    single_label(sample_info$response[sample_info$patient == p]), character(1))

  conflicted <- inconsistent_patients(comp, "response")
  cat("Patients with conflicting response labels across biopsies:",
      paste(conflicted, collapse = ", "), "\n")
  summary_rows <- list()

  ## 1. Composition, per patient and per biopsy
  patient_comp <- aggregate_to_patient(comp, cell_types)
  cat(sprintf("Patient-level composition: %d patients, %d excluded as %s\n",
              nrow(patient_comp), sum(patient_comp$response == CONFLICTING), CONFLICTING))
  groups <- list(
    per_patient = patient_comp,
    all_samples = comp,
    baseline = comp[comp$timepoint == "Pre", , drop = FALSE]
  )
  for (name in names(groups)) {
    table_out <- compare_fractions(groups[[name]], cell_types, "response", NR, R_)
    write.csv(table_out, file.path(tables_dir, sprintf("response_composition_%s.csv", name)),
              row.names = FALSE)
    n <- attr(table_out, "n")
    cat(sprintf("%s: NR n=%d, R n=%d, %d cell types at q < 0.1\n", name, n[[NR]], n[[R_]],
                sum(table_out$qvalue < 0.1, na.rm = TRUE)))
  }

  ## 2. TCF7+ CD8 T cells, per biopsy and per patient
  cd8 <- cd8_mask(obj)
  cat(sprintf("CD8 T cells: %d\n", sum(cd8)))
  tcf7 <- fraction_positive(m, "TCF7", cd8, meta$sample, a$min_cells_per_sample)
  write.csv(data.frame(sample = names(tcf7), tcf7_pos_fraction_of_cd8 = as.numeric(tcf7),
                       sample_info[names(tcf7), c("response", "patient", "timepoint")]),
            file.path(tables_dir, "tcf7_cd8_fraction_by_sample.csv"), row.names = FALSE)

  tcf7_pt <- tapply(as.numeric(tcf7), patients[names(tcf7)], mean)
  tcf7_pt <- tcf7_pt[patient_labels[names(tcf7_pt)] != CONFLICTING]
  units <- list(list(unit = "patient", values = tcf7_pt, labels = patient_labels),
                list(unit = "biopsy", values = setNames(as.numeric(tcf7), names(tcf7)),
                     labels = labels))
  for (u in units) {
    value_auc <- auc_score(u$values, u$labels, R_)
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      signature = "TCF7+ fraction of CD8 T cells (from the paper)", unit = u$unit,
      n = length(u$values), auc = value_auc, cross_validated = "not needed (no fitting)")
    png(file.path(figures_dir, sprintf("response_tcf7_cd8_per_%s.png", u$unit)),
        width = 700, height = 800, res = 150)
    g <- factor(u$labels[names(u$values)], levels = c(NR, R_))
    stripchart(as.numeric(u$values) ~ g, vertical = TRUE, method = "jitter", pch = 16,
               col = c("#2a78d6", "#eb6834"), ylab = "TCF7+ fraction of CD8 T cells",
               main = sprintf("TCF7+ CD8 T cells, per %s (AUC=%.2f)", u$unit, value_auc),
               cex.main = 0.9)
    dev.off()
  }

  ## 3 & 4. Differential expression and learned signatures
  for (name in c("all_immune", "cd8")) {
    mask <- if (name == "cd8") cd8 else rep(TRUE, ncol(m))
    title <- if (name == "cd8") "CD8 T cells" else "All immune cells"
    pbm <- pseudobulk_matrices(m[, mask, drop = FALSE], meta$sample[mask],
                               min_cells = a$min_cells_per_sample)
    rate <- detection_rate(pbm$detect, pbm$n_cells)
    kept_genes <- names(rate)[!is.na(rate) & rate >= a$min_frac_expressed]
    pb <- pbm$means[, kept_genes, drop = FALSE]
    lab <- labels[rownames(pb)]
    de <- pseudobulk_de(pb, lab, NR, R_)
    write.csv(data.frame(gene = rownames(de), de),
              file.path(tables_dir, sprintf("response_de_%s.csv", name)), row.names = FALSE)
    cat(sprintf("%s: %d samples, %d genes, %d genes at q < 0.1\n", name, nrow(pb), ncol(pb),
                sum(de$qvalue < 0.1, na.rm = TRUE)))

    png(file.path(figures_dir, sprintf("response_de_%s.png", name)),
        width = 1000, height = 850, res = 150)
    sig <- de$qvalue < 0.1
    plot(de$log2fc, -log10(de$pvalue), pch = 16, cex = 0.35,
         col = ifelse(!sig, "#e4e3dd", ifelse(de$log2fc > 0, "#eb6834", "#2a78d6")),
         xlab = sprintf("log2 fold change (%s vs %s)", R_, NR), ylab = "-log10 p-value",
         main = sprintf("Responder vs non-responder, %s (per-biopsy means)", title),
         cex.main = 0.9)
    abline(v = 0, col = "#6b6a64", lwd = 0.6)
    top <- head(rownames(de), 12)
    text(de[top, "log2fc"], -log10(de[top, "pvalue"]), top, cex = 0.5, pos = 4)
    dev.off()

    folds <- list(list(tag = "all", note = "all biopsies", idx = rownames(pbm$means)))
    if (length(conflicted)) {
      keep <- rownames(pbm$means)[!(patients[rownames(pbm$means)] %in% conflicted)]
      folds[[2]] <- list(tag = "no_conflict",
                         note = sprintf("excluding %d label-conflicted patients",
                                        length(conflicted)), idx = keep)
    }
    for (f in folds) {
      scores <- leave_one_patient_out_scores(
        pbm$means[f$idx, , drop = FALSE], labels[f$idx], patients[f$idx], NR, R_,
        a$signature_size, detect = pbm$detect, n_cells = pbm$n_cells,
        min_frac_expressed = a$min_frac_expressed)
      cv_auc <- auc_score(scores[!is.na(scores)], labels, R_)
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        signature = sprintf("Learned %d+%d gene signature, %s", a$signature_size,
                            a$signature_size, title),
        unit = sprintf("biopsy, %s", f$note), n = sum(!is.na(scores)), auc = cv_auc,
        cross_validated = "leave-one-patient-out")
      write.csv(data.frame(sample = names(scores), signature_score = as.numeric(scores),
                           sample_info[names(scores), c("response", "patient", "timepoint")]),
                file.path(tables_dir,
                          sprintf("response_signature_scores_%s_%s.csv", name, f$tag)),
                row.names = FALSE)
      png(file.path(figures_dir, sprintf("response_signature_%s_%s.png", name, f$tag)),
          width = 700, height = 800, res = 150)
      ok <- !is.na(scores)
      g <- factor(labels[names(scores)[ok]], levels = c(NR, R_))
      stripchart(as.numeric(scores[ok]) ~ g, vertical = TRUE, method = "jitter", pch = 16,
                 col = c("#2a78d6", "#eb6834"), ylab = "Signature score (held-out)",
                 main = sprintf("Learned signature, %s, %s\n(cross-validated AUC=%.2f)",
                                title, f$note, cv_auc), cex.main = 0.75)
      dev.off()
    }
  }

  summary_df <- do.call(rbind, summary_rows)
  write.csv(summary_df, file.path(tables_dir, "response_signature_summary.csv"),
            row.names = FALSE)
  print(summary_df, digits = 3)
}

main()
