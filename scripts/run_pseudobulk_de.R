#!/usr/bin/env Rscript
# CM pseudobulk DE model (Plan Step 7). The heavy CM aggregation over ~1M cells is done in
# Python (pseudobulk_cm.py, backed read) which writes a small genes x patients table; this script just
# runs limma-voom (design ~ 0 + disease + study + anatomy + sex, duplicateCorrelation on individual) on
# that small table -> no big in-memory object here.
#
# Usage:
#   run_pseudobulk_de.R --pb_counts <cm_pseudobulk_counts.tsv> --pb_meta <cm_pseudobulk_meta.tsv> \
#                       --out <de_dir>

suppressPackageStartupMessages({
  library(optparse); library(edgeR); library(limma)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--pb_counts"),
  make_option("--pb_meta"),
  make_option("--out"),
  make_option("--patient_col",  default = "individual"),
  make_option("--disease_col",  default = "disease"),
  make_option("--anatomy_col",  default = "anatomy"),
  make_option("--sex_col",      default = "sex"),
  make_option("--study_col",    default = "study")
)))
stopifnot(!is.null(opt$pb_counts), !is.null(opt$pb_meta), !is.null(opt$out))
dir.create(opt$out, showWarnings = FALSE, recursive = TRUE)

# ---- read the small pseudobulk tables (genes x patients + patient meta) ----
pb <- as.matrix(read.delim(opt$pb_counts, row.names = 1, check.names = FALSE))
pmeta <- read.delim(opt$pb_meta, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
pmeta <- pmeta[colnames(pb), , drop = FALSE]
message(sprintf("pseudobulk: %d genes x %d patients", nrow(pb), ncol(pb)))

# ---- limma-voom + duplicateCorrelation ----
dge <- DGEList(pb)
keep <- filterByExpr(dge, group = pmeta[[opt$disease_col]])
dge  <- dge[keep, , keep.lib.sizes = FALSE]
dge  <- calcNormFactors(dge)

# ---- metadata completeness guard ----
# Missing/NA disease or study silently corrupts the design (mis-assigned contrasts / hidden batch).
# Hard-fail on those; loudly warn on anatomy/sex (missing covariates become an explicit level, which
# is modeled, not silently dropped).
.na_like <- function(v) is.na(v) | trimws(as.character(v)) %in% c("", "NA", "nan", "None", "NaN")
for (mc in c(opt$disease_col, opt$study_col)) {
  bad <- .na_like(pmeta[[mc]])
  if (any(bad)) stop(sprintf("[DE] %d/%d patients have missing/NA '%s': refusing to build a corrupted "
                             , sum(bad), nrow(pmeta), mc),
                     "design. Complete cm_pseudobulk_meta.tsv (upstream: sample_meta.tsv) before DE.")
}
for (mc in c(opt$anatomy_col, opt$sex_col)) {
  bad <- .na_like(pmeta[[mc]])
  if (any(bad)) message(sprintf("[DE] WARNING: %d/%d patients have missing '%s'; modeled as an explicit "
                                , sum(bad), nrow(pmeta), mc), "level (verify this is intended).")
}

disease <- factor(make.names(as.character(pmeta[[opt$disease_col]])))
study   <- factor(make.names(as.character(pmeta[[opt$study_col]])))
anatomy <- factor(make.names(as.character(pmeta[[opt$anatomy_col]])))
sex     <- factor(make.names(as.character(pmeta[[opt$sex_col]])))
# Build the design from disease + only the covariates that actually vary (>=2 levels). A
# single-level covariate (e.g. one study in a within-cohort run, or an all-missing sex column)
# contributes no design columns and makes model.matrix() abort with "contrasts can be applied only
# to factors with 2 or more levels"; dropping it is a no-op for the fit. disease is always retained
# (the contrasts of interest are built from its coefficients).
covars  <- c("study", "anatomy", "sex")
kept    <- covars[vapply(covars, function(v) nlevels(get(v)) >= 2, logical(1))]
dropped <- setdiff(covars, kept)
if (length(dropped))
  message(sprintf("[DE] dropped single-level covariate(s) from design (uninformative): %s",
                  paste(dropped, collapse = ", ")))
design  <- model.matrix(as.formula(paste("~", paste(c("0", "disease", kept), collapse = " + "))))

# ---- rank / estimability guard ----
# If a disease is present in only one study, its disease coefficient is perfectly aliased with that
# study coefficient (e.g. CS == neyazi). limma would silently NA-drop an aliased coefficient and
# emit normal-looking topTables, so the resulting disease contrasts are batch-confounded, NOT
# batch-clean. Detect this explicitly and refuse to run silently.
ne <- limma::nonEstimable(design)
design_status <- "full_rank"                 # written into every output file for provenance
single_study <- character(0)
if (!limma::is.fullrank(design) || !is.null(ne)) {
  design_status <- "rank_deficient_EXPLORATORY_confounded"
  ct <- table(as.character(disease), as.character(study))
  single_study <- rownames(ct)[rowSums(ct > 0) == 1]
  message("[DE] *** RANK-DEFICIENT DESIGN ***  rank=", qr(design)$rank, " of ", ncol(design), " columns")
  if (!is.null(ne)) message("[DE] non-estimable coefficients: ", paste(ne, collapse = ", "))
  if (length(single_study))
    message("[DE] single-study diseases (aliased with their study, contrasts CONFOUNDED): ",
            paste(single_study, collapse = ", "))
  message("[DE] These disease contrasts are EXPLORATORY, not batch-controlled; do not present as ",
          "specificity evidence. Batch-clean specificity must come from within-study/independent data.")
  if (Sys.getenv("ALLOW_RANK_DEFICIENT", "0") != "1")
    stop("[DE] design not full rank. Set ALLOW_RANK_DEFICIENT=1 to proceed on confounded (exploratory) ",
         "contrasts; the output must then be labeled exploratory.")
}

v <- voom(dge, design)
# duplicateCorrelation only makes sense when some individual has >1 sample. With one pseudobulk
# sample per patient every block is unique -> dupCor returns NaN and breaks lmFit; fall back to
# plain lmFit in that case (also covers real cohorts where each donor contributes a single sample).
block <- colnames(pb)
if (any(duplicated(block))) {
  corfit <- duplicateCorrelation(v, design, block = block)
  fit <- lmFit(v, design, block = block, correlation = corfit$consensus)
} else {
  message("[DE] no repeated individuals -> plain lmFit (duplicateCorrelation skipped)")
  fit <- lmFit(v, design)
}

# ---- procurement-matched primary contrasts ----
lvl <- colnames(design)
mk  <- function(a, b) {
  ca <- paste0("disease", make.names(a)); cb <- paste0("disease", make.names(b))
  if (!(ca %in% lvl) || !(cb %in% lvl)) return(NULL)
  setNames((lvl == ca) - (lvl == cb), NULL)
}
contrasts <- list(
  CS_vs_DCM          = mk("CS", "DCM"),
  CS_vs_ICM_Simonson = mk("CS", "ICM"),
  CS_vs_ARVC         = mk("CS", "ARVC"),
  CS_vs_HCM          = mk("CS", "HCM"),
  CS_vs_NF           = mk("CS", "NF")
)
contrasts <- contrasts[!sapply(contrasts, is.null)]

for (nm in names(contrasts)) {
  cm_fit <- contrasts.fit(fit, contrasts[[nm]])
  cm_fit <- eBayes(cm_fit)
  tt <- topTable(cm_fit, number = Inf, sort.by = "P")
  # Provenance stamp IN THE FILE: anyone reading the topTable months later must see the
  # design status without consulting stderr. If a CS contrast is confounded, mark it per-row.
  contrast_confounded <- design_status != "full_rank" &&
    any(sapply(single_study, function(d) grepl(d, nm, ignore.case = TRUE)))
  tt$design_status <- design_status
  tt$contrast_confounded <- contrast_confounded
  out <- file.path(opt$out, sprintf("DE_%s.tsv", nm))
  write.table(tt, out, sep = "\t", quote = FALSE, col.names = NA)
  message(sprintf("wrote %s (%d genes, %d FDR<0.05; design=%s, contrast_confounded=%s)",
                  out, nrow(tt), sum(tt$adj.P.Val < 0.05, na.rm = TRUE),
                  design_status, contrast_confounded))
}
message("Pseudobulk DE complete. NOTE: combine with pydeseq2/apeglm shrinkage + >=75% per-study ",
        "concordance + procurement gate for the final robust-DEG list (Plan Step 7).")
