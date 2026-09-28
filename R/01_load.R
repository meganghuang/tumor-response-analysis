## Load GSE120575 from the same raw GEO text files the Python path reads.
##
## The expression file is genes x cells log2(TPM+1) with two header lines (cell titles, then
## sample IDs) and a trailing tab on every data row. Values are converted to ln(TPM+1), the
## scale Seurat's `data` layer is expected to hold. Rows are read in blocks and sparsified
## immediately: the full matrix is ~55,000 x 16,291, which would need about 7 GB dense.

source(file.path("R", "00_common.R"))
suppressPackageStartupMessages(library(data.table))

EXPR_FILE <- "GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz"
META_FILE <- "GSE120575_patient_ID_single_cells.txt.gz"

#' Per-cell metadata indexed by cell title. The GEO file is a submission template: comment
#' lines, then one row per cell, then trailing template sections, so rows are selected by
#' the "Sample N" pattern rather than by position.
read_metadata <- function(path) {
  con <- gzfile(path, "rt")
  lines <- readLines(con, warn = FALSE)
  close(con)
  header_idx <- grep("^Sample name\t", lines)[1]
  stopifnot(!is.na(header_idx))
  header <- trimws(strsplit(lines[header_idx], "\t")[[1]])
  body <- lines[grepl("^Sample [0-9]+\t", lines)]
  fields <- strsplit(body, "\t")
  col <- function(pattern) {
    j <- grep(pattern, header)[1]
    stopifnot(!is.na(j))
    trimws(vapply(fields, function(f) if (length(f) >= j) f[j] else NA_character_, character(1)))
  }
  meta <- data.frame(
    sample = col("^characteristics: patinet ID"),
    response = col("^characteristics: response"),
    therapy = col("^characteristics: therapy"),
    row.names = col("^title$"),
    stringsAsFactors = FALSE
  )
  parts <- strsplit(meta$sample, "_")
  meta$timepoint <- vapply(parts, `[`, character(1), 1L)
  meta$patient <- vapply(parts, `[`, character(1), 2L)
  meta
}

#' Decompress once to a plain text file next to the archive, and return that path.
#' fread() can only read .gz through R.utils, and reading blocks from a plain file avoids
#' re-inflating the whole archive on every block.
ensure_uncompressed <- function(path) {
  out <- sub("\\.gz$", "", path)
  if (file.exists(out) && file.size(out) > file.size(path)) return(out)
  cat("Decompressing", basename(path), "\n"); flush.console()
  gz <- gzfile(path, "rb")
  con <- file(out, "wb")
  repeat {
    chunk <- readBin(gz, "raw", 8e6)
    if (!length(chunk)) break
    writeBin(chunk, con)
  }
  close(gz); close(con)
  out
}

#' Genes x cells sparse matrix of ln(TPM+1).
read_expression <- function(path, block = 4000L) {
  path <- ensure_uncompressed(path)
  header <- strsplit(readLines(path, 1L), "\t")[[1]]
  cells <- trimws(header[-1])
  cells <- cells[nzchar(cells)]
  n_cells <- length(cells)

  blocks <- list()
  genes <- character(0)
  skip <- 2L
  repeat {
    dt <- data.table::fread(path, sep = "\t", skip = skip, nrows = block, header = FALSE,
                            showProgress = FALSE)
    if (nrow(dt) == 0L) break
    stopifnot(ncol(dt) >= n_cells + 1L)
    genes <- c(genes, as.character(dt[[1L]]))
    values <- as.matrix(dt[, 2:(n_cells + 1L), with = FALSE])
    values[is.na(values)] <- 0
    ## log2(TPM+1) * ln(2) = ln(TPM+1)
    blocks[[length(blocks) + 1L]] <- Matrix::Matrix(values * log(2), sparse = TRUE)
    skip <- skip + nrow(dt)
    cat(sprintf("  read %d genes\n", length(genes))); flush.console()
    if (nrow(dt) < block) break
  }
  m <- do.call(rbind, blocks)
  rownames(m) <- make.unique(genes)
  colnames(m) <- cells
  m
}

load_sade_feldman <- function(raw_dir) {
  m <- read_expression(file.path(raw_dir, EXPR_FILE))
  meta <- read_metadata(file.path(raw_dir, META_FILE))
  if (anyDuplicated(colnames(m))) stop("duplicate cell titles in the expression matrix")
  if (anyDuplicated(rownames(meta))) stop("duplicate cell titles in the metadata")
  missing <- setdiff(colnames(m), rownames(meta))
  if (length(missing)) {
    stop(sprintf("%d cells have no metadata, e.g. %s", length(missing),
                 paste(head(missing, 3), collapse = ", ")))
  }
  list(matrix = m, meta = meta[colnames(m), , drop = FALSE])
}

if (sys.nframe() == 0L) {
  cfg <- load_cfg()
  raw_dir <- project_path(cfg$data$raw_dir)
  processed_dir <- project_path(cfg$data$processed_dir)
  cache <- file.path(processed_dir, "00_raw_seurat.rds")
  if (file.exists(cache)) {
    cat("Already parsed:", cache, "\n")
  } else {
    data <- load_sade_feldman(raw_dir)
    cat(sprintf("Loaded: %d cells x %d genes, %d samples\n", ncol(data$matrix),
                nrow(data$matrix), length(unique(data$meta$sample))))
    saveRDS(data, cache, compress = FALSE)
  }
}
