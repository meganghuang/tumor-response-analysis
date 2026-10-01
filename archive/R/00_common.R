## Shared configuration, paths and statistics for the Seurat path.
##
## The statistics mirror src/tumor_response/: per-sample or per-patient values, Mann-Whitney U
## (unpaired) or Wilcoxon signed-rank (paired), Benjamini-Hochberg q-values, and
## leave-one-patient-out cross-validation with everything refitted inside the fold.
##
## Run from the project root, or set TUMOR_RESPONSE_ROOT.

suppressPackageStartupMessages({
  library(Matrix)
  library(yaml)
})

project_root <- function() {
  normalizePath(Sys.getenv("TUMOR_RESPONSE_ROOT", unset = getwd()), mustWork = TRUE)
}

load_cfg <- function() yaml::read_yaml(file.path(project_root(), "config", "config.yaml"))

project_path <- function(relative) {
  path <- file.path(project_root(), relative)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, mustWork = TRUE)
}

## Output directories for the R path, kept separate so neither path overwrites the other.
r_figures_dir <- function(cfg) project_path(file.path(cfg$output$figures_dir, "seurat"))
r_tables_dir <- function(cfg) project_path(file.path(cfg$output$tables_dir, "seurat"))

CONFLICTING <- "Conflicting"

## --- statistics ---------------------------------------------------------------------------

#' The one label shared by a patient's samples, or "Conflicting" if they disagree.
single_label <- function(values) {
  u <- unique(as.character(values[!is.na(values)]))
  if (length(u) == 1L) u else CONFLICTING
}

#' P(a random value from `b` exceeds a random value from `a`), i.e. the Mann-Whitney U
#' statistic divided by n_a * n_b, which is the ROC AUC.
auc_u <- function(b, a) {
  b <- b[!is.na(b)]; a <- a[!is.na(a)]
  if (!length(a) || !length(b)) return(NA_real_)
  as.numeric(suppressWarnings(wilcox.test(b, a, exact = FALSE))$statistic) /
    (length(a) * length(b))
}

#' AUC of a score against a two-level label, for the named positive class.
auc_score <- function(scores, labels, positive) {
  labels <- labels[names(scores)]
  ok <- !is.na(scores) & !is.na(labels)
  auc_u(scores[ok & labels == positive], scores[ok & labels != positive])
}

#' Unpaired comparison of one value per row between two groups, for several columns.
compare_fractions <- function(df, cell_types, group_col, group_a, group_b) {
  a <- df[df[[group_col]] == group_a, , drop = FALSE]
  b <- df[df[[group_col]] == group_b, , drop = FALSE]
  out <- do.call(rbind, lapply(cell_types, function(ct) {
    x <- a[[ct]]; y <- b[[ct]]
    test <- suppressWarnings(wilcox.test(y, x, exact = FALSE))
    data.frame(cell_type = ct, mean_a = mean(x), mean_b = mean(y), diff = mean(y) - mean(x),
               auc = as.numeric(test$statistic) / (length(x) * length(y)),
               pvalue = test$p.value)
  }))
  names(out)[2:3] <- paste0("mean_", c(group_a, group_b))
  out$qvalue <- p.adjust(out$pvalue, method = "BH")
  attr(out, "n") <- setNames(c(nrow(a), nrow(b)), c(group_a, group_b))
  out[order(out$pvalue), ]
}

#' Paired before/after comparison for patients sampled at both timepoints. Patients with
#' several lesions at one timepoint are averaged first.
paired_fractions <- function(comp, cell_types, patient_col = "patient",
                             time_col = "timepoint", before = "Pre", after = "Post") {
  per <- aggregate(comp[cell_types],
                   by = list(patient = comp[[patient_col]], timepoint = comp[[time_col]]),
                   FUN = mean)
  pre <- per[per$timepoint == before, ]
  post <- per[per$timepoint == after, ]
  shared <- sort(intersect(pre$patient, post$patient))
  pre <- pre[match(shared, pre$patient), ]
  post <- post[match(shared, post$patient), ]
  out <- do.call(rbind, lapply(cell_types, function(ct) {
    d <- post[[ct]] - pre[[ct]]
    test <- if (any(d != 0)) suppressWarnings(
      wilcox.test(post[[ct]], pre[[ct]], paired = TRUE, exact = FALSE)) else NULL
    data.frame(cell_type = ct, mean_Pre = mean(pre[[ct]]), mean_Post = mean(post[[ct]]),
               mean_change = mean(d), n_increase = sum(d > 0), n_decrease = sum(d < 0),
               pvalue = if (is.null(test)) NA_real_ else test$p.value)
  }))
  out$qvalue <- p.adjust(out$pvalue, method = "BH")
  attr(out, "n_patients") <- length(shared)
  attr(out, "pre") <- pre
  attr(out, "post") <- post
  out[order(out$pvalue), ]
}

#' Fraction of each cell type per sample, plus that sample's annotations.
#' Proportions are per sample, so samples with more cells do not dominate.
composition_by_sample <- function(meta, sample_col = "sample", cell_type_col = "cell_type",
                                  annotations = c("patient", "timepoint", "response", "therapy")) {
  counts <- table(meta[[sample_col]], meta[[cell_type_col]])
  fractions <- as.data.frame.matrix(counts / rowSums(counts))
  samples <- rownames(fractions)
  for (col in annotations) {
    fractions[[col]] <- vapply(samples, function(s)
      single_label(meta[meta[[sample_col]] == s, col]), character(1))
  }
  fractions[[sample_col]] <- samples
  fractions
}

#' Per sample, the fraction of cells in `mask` with the gene detected (e.g. TCF7+ of CD8).
fraction_positive <- function(m, gene, mask, samples, min_cells = 10) {
  stopifnot(gene %in% rownames(m))
  x <- m[gene, mask, drop = TRUE] > 0
  s <- as.character(samples[mask])
  frac <- tapply(x, s, mean)
  sizes <- tapply(x, s, length)
  frac[sizes >= min_cells]
}

#' Samples from patients biopsied at only one timepoint, so Pre and Post share no patient.
single_timepoint_samples <- function(comp, patient_col = "patient", time_col = "timepoint") {
  n <- tapply(comp[[time_col]], comp[[patient_col]], function(x) length(unique(x)))
  comp[comp[[patient_col]] %in% names(n)[n == 1], , drop = FALSE]
}

#' One row per patient: mean cell-type fraction, with conflicting labels marked.
aggregate_to_patient <- function(comp, cell_types, patient_col = "patient",
                                 annotations = c("response", "therapy")) {
  patients <- sort(unique(comp[[patient_col]]))
  out <- as.data.frame(do.call(rbind, lapply(patients, function(p) {
    colMeans(comp[comp[[patient_col]] == p, cell_types, drop = FALSE])
  })))
  rownames(out) <- patients
  for (col in intersect(annotations, names(comp))) {
    out[[col]] <- vapply(patients, function(p)
      single_label(comp[comp[[patient_col]] == p, col]), character(1))
  }
  out[[patient_col]] <- patients
  out
}

#' Patients whose biopsies carry more than one distinct value of `col`.
inconsistent_patients <- function(comp, col, patient_col = "patient") {
  n <- tapply(comp[[col]], comp[[patient_col]], function(x) length(unique(x)))
  sort(names(n)[n > 1])
}

#' Per-gene Mann-Whitney U test between samples of two groups, on a samples x genes matrix.
pseudobulk_de <- function(pb, labels, group_a, group_b) {
  a <- pb[labels == group_a, , drop = FALSE]
  b <- pb[labels == group_b, , drop = FALSE]
  na <- nrow(a); nb <- nrow(b)
  ranks <- apply(rbind(a, b), 2, rank)
  u_b <- colSums(ranks[(na + 1):(na + nb), , drop = FALSE]) - nb * (nb + 1) / 2
  pvalue <- vapply(seq_len(ncol(pb)), function(j) suppressWarnings(
    wilcox.test(b[, j], a[, j], exact = FALSE)$p.value), numeric(1))
  mean_a <- colMeans(a); mean_b <- colMeans(b)
  out <- data.frame(mean_a = mean_a, mean_b = mean_b,
                    log2fc = (mean_b - mean_a) / log(2), auc = u_b / (na * nb),
                    pvalue = pvalue, row.names = colnames(pb))
  names(out)[1:2] <- paste0("mean_", c(group_a, group_b))
  out$qvalue <- p.adjust(out$pvalue, method = "BH")
  out[order(out$pvalue), ]
}

#' The n most significant genes up and down in group_b. Fewer if fewer genes move that way.
top_genes <- function(de, n) {
  list(up = head(rownames(de)[de$log2fc > 0], n),
       down = head(rownames(de)[de$log2fc < 0], n))
}

#' Mean z-score of `up` minus mean z-score of `down`, z-scored against `ref`.
signature_score <- function(pb, up, down, ref) {
  mu <- colMeans(ref)
  sdev <- apply(ref, 2, sd)
  sdev[sdev == 0] <- NA_real_
  z <- sweep(sweep(pb, 2, mu, "-"), 2, sdev, "/")
  score <- if (length(up)) rowMeans(z[, up, drop = FALSE]) else 0
  if (length(down)) score <- score - rowMeans(z[, down, drop = FALSE])
  score
}

#' Fraction of cells with each gene detected, pooled over the named samples.
detection_rate <- function(detect, n_cells, samples = rownames(detect)) {
  d <- detect[samples, , drop = FALSE]
  w <- n_cells[samples]
  colSums(d * w) / sum(w)
}

#' Score every sample with a signature learned without any sample from its patient.
#' The detection filter, the gene selection and the z-reference are fitted on training
#' samples only, so the held-out patient never influences its own score.
leave_one_patient_out_scores <- function(pb, labels, patients, group_a, group_b, n_genes,
                                         detect = NULL, n_cells = NULL,
                                         min_frac_expressed = 0) {
  scores <- setNames(rep(NA_real_, nrow(pb)), rownames(pb))
  for (patient in unique(patients)) {
    test <- patients == patient
    train <- !test
    if (length(unique(labels[train])) < 2L) next
    sub <- pb
    if (!is.null(detect) && !is.null(n_cells) && min_frac_expressed > 0) {
      rate <- detection_rate(detect, n_cells, rownames(pb)[train])
      sub <- pb[, names(rate)[!is.na(rate) & rate >= min_frac_expressed], drop = FALSE]
    }
    de <- pseudobulk_de(sub[train, , drop = FALSE], labels[train], group_a, group_b)
    g <- top_genes(de, n_genes)
    scores[test] <- signature_score(sub[test, , drop = FALSE], g$up, g$down,
                                    ref = sub[train, , drop = FALSE])
  }
  scores
}

#' One-vs-rest Wilcoxon rank-sum marker test over every gene, vectorised.
#'
#' Same statistic as Seurat's `FindAllMarkers(test.use = "wilcox")` and scanpy's
#' `rank_genes_groups(method = "wilcoxon")`: for each group, each gene's expression is ranked
#' across all cells and the rank sum inside the group gives U, converted to a z-score with
#' the tie correction and then to a two-sided p-value by normal approximation.
#'
#' Seurat's own implementation is not usable at this size on this platform: without presto it
#' falls back to calling `wilcox.test` per gene per group, which does not finish in hours for
#' 14,505 cells. Ranking each gene once and reusing those ranks for every group turns the
#' same computation into minutes. Genes are processed in blocks to bound memory.
#'
#' @param m genes x cells matrix (the `data` layer: ln(TPM+1)).
#' @param groups per-cell group label.
#' @param block number of genes to rank at a time.
#' @return data.frame with one row per gene per group, ordered by z descending within group.
rank_sum_markers <- function(m, groups, block = 2000L, min_pct = 0.1) {
  groups <- as.character(groups)
  stopifnot(length(groups) == ncol(m))
  levels_g <- sort(unique(groups))
  n <- ncol(m)
  indicator <- vapply(levels_g, function(g) as.numeric(groups == g), numeric(n))
  n_in <- colSums(indicator)

  out <- vector("list", 0L)
  starts <- seq(1L, nrow(m), by = block)
  for (s in starts) {
    idx <- s:min(s + block - 1L, nrow(m))
    dense <- as.matrix(m[idx, , drop = FALSE])
    ## rank() and a run-length pass over the sorted values give the average ranks and the
    ## tie correction in one sweep. table() per gene would be correct but is ~100x slower,
    ## which at 35,000 genes is the difference between minutes and hours.
    ranks <- t(apply(dense, 1, rank, ties.method = "average"))
    tie_term <- apply(dense, 1, function(x) {
      lens <- rle(sort.int(x, method = "quick"))$lengths
      sum(lens^3 - lens)
    })
    linear <- expm1(dense)
    detected <- dense > 0
    for (j in seq_along(levels_g)) {
      w <- indicator[, j]
      n1 <- n_in[j]; n2 <- n - n1
      rank_sum <- as.numeric(ranks %*% w)
      u <- rank_sum - n1 * (n1 + 1) / 2
      mu <- n1 * n2 / 2
      sigma <- sqrt((n1 * n2 / 12) * ((n + 1) - tie_term / (n * (n - 1))))
      z <- ifelse(sigma > 0, (u - mu) / sigma, 0)
      mean_in <- as.numeric(linear %*% w) / n1
      mean_out <- as.numeric(linear %*% (1 - w)) / n2
      pct_in <- as.numeric(detected %*% w) / n1
      pct_out <- as.numeric(detected %*% (1 - w)) / n2
      out[[length(out) + 1L]] <- data.frame(
        cluster = levels_g[j], gene = rownames(dense),
        avg_log2FC = log2((mean_in + 1) / (mean_out + 1)),
        pct.1 = pct_in, pct.2 = pct_out, z = z,
        p_val = 2 * pnorm(-abs(z)), stringsAsFactors = FALSE)
    }
    cat(sprintf("  markers: %d / %d genes\n", max(idx), nrow(m))); flush.console()
  }
  res <- do.call(rbind, out)
  res <- res[res$pct.1 >= min_pct | res$pct.2 >= min_pct, , drop = FALSE]
  res$p_val_adj <- ave(res$p_val, res$cluster, FUN = function(p) p.adjust(p, method = "BH"))
  res[order(res$cluster, -res$z), ]
}

#' Per-sample mean log expression, detection rate, and cell count over all genes.
#' `m` is genes x cells; `samples` is a per-cell sample label.
pseudobulk_matrices <- function(m, samples, min_cells = 10) {
  samples <- as.character(samples)
  stopifnot(length(samples) == ncol(m))
  sizes <- table(samples)
  kept <- names(sizes)[sizes >= min_cells]
  ## Built explicitly rather than with sparse.model.matrix, which would silently drop the
  ## cells of excluded samples and misalign the rows against the columns of `m`.
  f <- factor(samples, levels = kept)
  use <- !is.na(f)
  onehot <- Matrix::sparseMatrix(i = which(use), j = as.integer(f[use]), x = 1,
                                 dims = c(length(samples), length(kept)),
                                 dimnames = list(NULL, kept))
  n_cells <- setNames(Matrix::colSums(onehot), kept)
  sums <- as.matrix(Matrix::t(onehot) %*% Matrix::t(m))
  detected <- as.matrix(Matrix::t(onehot) %*% Matrix::t(1 * (m > 0)))
  list(means = sums / n_cells, detect = detected / n_cells, n_cells = n_cells)
}
