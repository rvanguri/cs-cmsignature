#!/usr/bin/env Rscript
# Uniform DecontX pass (Plan Step 2b), run PER BATCH (per sample / per donor) and concatenated.
#   - GEO + Liu: one batch per manifest sample (CellBender-filtered h5, or raw filtered matrix).
#   - Reichart:  the CELLxGENE 'All cells' raw.X, split into donor batches (880k nuclei OOMs whole).
# Per-batch keeps peak memory bounded (whole-dataset DecontX builds a UMAP over ALL cells -> OOM on
# Reichart) and is DecontX's own recommended usage. Identical default params across every batch and
# dataset = the consistency fix. Output: integer decontaminated counts as .h5ad.
#
# Sample selection is MANIFEST-DRIVEN. Per-cell `sample_id` is written to obs (join key for sample_meta
# in run_qc). Reichart additionally carries disease/anatomy/sex/individual from its own CELLxGENE obs.
#
# Usage:
#   run_decontx.R --dataset <name> --manifest <manifest.tsv> --cellbender_dir <dir> \
#                 --reichart_dir <dir> --out <out.h5ad>

if (Sys.getenv("RETICULATE_PYTHON") == "" && Sys.getenv("CONDA_PREFIX") != "")
  Sys.setenv(RETICULATE_PYTHON = file.path(Sys.getenv("CONDA_PREFIX"), "bin", "python"))

suppressPackageStartupMessages({
  library(optparse); library(celda); library(SingleCellExperiment)
  library(Matrix); library(anndata)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--dataset"), make_option("--manifest"),
  make_option("--cellbender_dir"), make_option("--reichart_dir"), make_option("--out")
)))
stopifnot(!is.null(opt$dataset), !is.null(opt$out))
message(sprintf("[DecontX] dataset=%s", opt$dataset))

# light normalizers for Reichart CELLxGENE metadata -> shared vocab
map_region <- function(x) {
  x <- tolower(as.character(x))
  ifelse(grepl("left ventricle", x), "LV",
  ifelse(grepl("right ventricle", x), "RV",
  ifelse(grepl("interventricular|septum", x), "IVS",
  ifelse(grepl("apex", x), "apex", toupper(x)))))
}
map_sex <- function(x) {
  x <- tolower(as.character(x)); ifelse(x == "male", "M", ifelse(x == "female", "F", "NA"))
}

# ---- build a NAMED LIST of per-batch (genes x cells) matrices; Reichart also returns obs metadata ----
REICHART_OBS <- NULL   # filled for reichart; a data.frame aligned to the reichart matrix columns
load_batches <- function(ds, manifest_path, cb_dir, reichart_dir) {
  if (ds == "reichart") {
    ad  <- anndata::read_h5ad(file.path(reichart_dir, "reichart_cxg.h5ad"))
    raw <- ad$raw
    # reticulate returns X WITHOUT dimnames -> gene names live in var_names, must be attached
    if (!is.null(raw)) {
      mat   <- raw$X
      genes <- tryCatch(raw$var_names, error = function(e) NULL)
      if (is.null(genes)) genes <- rownames(raw$var)
    } else {
      mat   <- ad$X
      genes <- ad$var_names
    }
    cnt <- t(as(mat, "CsparseMatrix"))                    # genes x cells
    if (length(genes) != nrow(cnt))
      stop(sprintf("reichart gene-name length %d != matrix rows %d", length(genes), nrow(cnt)))
    rownames(cnt) <- genes                                # <-- the fix: real gene ids on the rows
    colnames(cnt) <- paste0("reichart|", seq_len(ncol(cnt)))
    obs <- ad$obs
    # stash normalized metadata for the whole dataset (indexed by the same column order)
    dcol <- intersect(c("disease","disease_state"), colnames(obs))
    rcol <- intersect(c("tissue","anatomy","region"), colnames(obs))
    scol <- intersect(c("sex"), colnames(obs))
    icol <- intersect(c("donor_id","donor","individual"), colnames(obs))
    ccol <- intersect(c("cell_type","cell_type_ontology","celltype","author_cell_type"), colnames(obs))
    REICHART_OBS <<- data.frame(
      row.names = colnames(cnt),
      disease    = if (length(dcol)) as.character(obs[[dcol[1]]]) else "NA",
      anatomy    = if (length(rcol)) map_region(obs[[rcol[1]]])  else "NA",
      sex        = if (length(scol)) map_sex(obs[[scol[1]]])     else "NA",
      individual = if (length(icol)) as.character(obs[[icol[1]]]) else "reichart",
      cell_type  = if (length(ccol)) as.character(obs[[ccol[1]]]) else "Unknown",  # scANVI reference labels
      stringsAsFactors = FALSE)
    # batch by donor to keep DecontX tractable
    if (length(icol) == 0) {
      message("[DecontX] reichart: no donor column -> single batch (may OOM)")
      return(list(reichart = cnt))
    }
    b <- as.character(obs[[icol[1]]])
    out <- lapply(split(seq_len(ncol(cnt)), b), function(idx) cnt[, idx, drop = FALSE])
    message(sprintf("[DecontX] reichart: %d donor batches", length(out)))
    return(out)
  }

  man  <- read.delim(manifest_path, stringsAsFactors = FALSE)
  rows <- man[man$dataset == ds, , drop = FALSE]
  if (nrow(rows) == 0) stop(sprintf("no manifest rows for dataset '%s'", ds))
  message(sprintf("[DecontX] %s: %d samples from manifest", ds, nrow(rows)))
  out <- list()
  for (i in seq_len(nrow(rows))) {
    sid <- rows$sample_id[i]
    cb  <- file.path(cb_dir, paste0(sid, "_cellbender_filtered.h5"))
    src <- if (file.exists(cb)) cb else rows$raw_path[i]
    if (!file.exists(cb))
      message(sprintf("[DecontX] %s: no CellBender output -> filtered matrix %s", sid, src))
    m  <- DropletUtils::read10xCounts(src)
    cm <- as(counts(m), "CsparseMatrix")                  # coerce DelayedMatrix/dgCMatrix uniformly
    colnames(cm) <- paste(sid, colData(m)$Barcode, sep = "|")
    out[[sid]] <- cm
  }
  out
}

batches <- load_batches(opt$dataset, opt$manifest, opt$cellbender_dir, opt$reichart_dir)

# ---- DecontX per batch, dropping empty cells; collect decontaminated matrices ----
dec_list <- list(); contam_all <- numeric(0)
for (nm in names(batches)) {
  cm  <- batches[[nm]]; batches[[nm]] <- NULL       # free the input copy as we go
  lib <- Matrix::colSums(cm); keep <- lib > 0
  if (any(!keep)) cm <- cm[, keep, drop = FALSE]
  if (ncol(cm) < 2) { message(sprintf("[DecontX] batch %s <2 cells -> skip", nm)); next }
  gn <- rownames(cm); cn <- colnames(cm)            # decontXcounts drops dimnames -> keep them
  sce <- SingleCellExperiment(assays = list(counts = cm))
  set.seed(1)
  sce <- tryCatch(decontX(sce), error = function(e) {
           message(sprintf("[DecontX] batch %s failed (%s) -> passing raw counts through",
                           nm, conditionMessage(e))); NULL })
  if (is.null(sce)) { dec <- cm; ctm <- rep(NA_real_, ncol(cm)) }
  else {
    dec <- as(round(decontXcounts(sce)), "CsparseMatrix")   # keep SPARSE (dense cbind would OOM)
    ctm <- sce$decontX_contamination
  }
  rownames(dec) <- gn; colnames(dec) <- cn          # restore gene/cell names (fixes "no common genes")
  dec_list[[nm]] <- dec; contam_all <- c(contam_all, ctm)
  rm(sce, cm); gc()
}
if (length(dec_list) == 0) stop("no batches produced output")

# ---- concatenate on the common gene set ----
common  <- Reduce(intersect, lapply(dec_list, rownames))
if (length(common) == 0) stop("no common genes across batches")
dec_all <- do.call(cbind, lapply(dec_list, function(x) x[common, , drop = FALSE]))
cells   <- colnames(dec_all)

# ---- assemble obs: sample_id + contamination (+ Reichart's own metadata) ----
if (opt$dataset == "reichart") {
  meta <- REICHART_OBS[cells, , drop = FALSE]
  obs  <- data.frame(row.names = cells,
                     sample_id  = "reichart",
                     disease    = meta$disease, anatomy = meta$anatomy,
                     sex        = meta$sex,     individual = meta$individual,
                     cell_type  = meta$cell_type,      # scANVI reference labels (mapped in run_scanvi)
                     decontX_contamination = contam_all, stringsAsFactors = FALSE)
} else {
  obs <- data.frame(row.names = cells,
                    sample_id = sub("\\|.*$", "", cells),        # sample id is the prefix before '|'
                    decontX_contamination = contam_all, stringsAsFactors = FALSE)
}

ad_out <- AnnData(X = t(as(dec_all, "CsparseMatrix")),
                  obs = obs, var = data.frame(row.names = rownames(dec_all)))
ad_out$write_h5ad(opt$out)
message(sprintf("[DecontX] wrote %s (%d cells x %d genes, %d batches)",
                opt$out, ncol(dec_all), nrow(dec_all), length(dec_list)))
