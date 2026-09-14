#!/usr/bin/env Rscript
# Patient-blocked variant of run_pseudobulk_de.R (robustness check).
# Identical to the original EXCEPT the duplicateCorrelation block is derived at the PATIENT level
# (not the per-sample column name), so neyazi CS patients contributing multiple region/replicate
# pseudobulk columns are modeled as repeated measures instead of independent observations.
# Outputs to a SEPARATE dir so the committed de_cm/ tables are untouched.
#
# Usage:
#   run_pseudobulk_de_patientblock.R --pb_counts <cm_pseudobulk_counts.tsv> \
#       --pb_meta <cm_pseudobulk_meta.tsv> --out <de_dir_patientblock>

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

pb <- as.matrix(read.delim(opt$pb_counts, row.names = 1, check.names = FALSE))
pmeta <- read.delim(opt$pb_meta, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
pmeta <- pmeta[colnames(pb), , drop = FALSE]
message(sprintf("pseudobulk: %d genes x %d samples", nrow(pb), ncol(pb)))

dge <- DGEList(pb)
keep <- filterByExpr(dge, group = pmeta[[opt$disease_col]])
dge  <- dge[keep, , keep.lib.sizes = FALSE]
dge  <- calcNormFactors(dge)

.na_like <- function(v) is.na(v) | trimws(as.character(v)) %in% c("", "NA", "nan", "None", "NaN")
for (mc in c(opt$disease_col, opt$study_col)) {
  bad <- .na_like(pmeta[[mc]])
  if (any(bad)) stop(sprintf("[DE] %d/%d samples have missing/NA '%s'", sum(bad), nrow(pmeta), mc))
}

disease <- factor(make.names(as.character(pmeta[[opt$disease_col]])))
study   <- factor(make.names(as.character(pmeta[[opt$study_col]])))
anatomy <- factor(make.names(as.character(pmeta[[opt$anatomy_col]])))
sex     <- factor(make.names(as.character(pmeta[[opt$sex_col]])))
covars  <- c("study", "anatomy", "sex")
kept    <- covars[vapply(covars, function(v) nlevels(get(v)) >= 2, logical(1))]
design  <- model.matrix(as.formula(paste("~", paste(c("0", "disease", kept), collapse = " + "))))

design_status <- "full_rank"
single_study <- character(0)
ne <- limma::nonEstimable(design)
if (!limma::is.fullrank(design) || !is.null(ne)) {
  design_status <- "rank_deficient_EXPLORATORY_confounded"
  ct <- table(as.character(disease), as.character(study))
  single_study <- rownames(ct)[rowSums(ct > 0) == 1]
  message("[DE] *** RANK-DEFICIENT DESIGN ***  rank=", qr(design)$rank, " of ", ncol(design))
  if (length(single_study))
    message("[DE] single-study diseases (confounded): ", paste(single_study, collapse = ", "))
  if (Sys.getenv("ALLOW_RANK_DEFICIENT", "0") != "1")
    stop("[DE] design not full rank. Set ALLOW_RANK_DEFICIENT=1 to proceed (exploratory).")
}

# ---- PATIENT-LEVEL block (the fix) ----
# neyazi column names are GSM<acc>_<patient>_<region>_<rep>; patient = 2nd underscore token.
# every other study contributes exactly one pseudobulk column per donor -> patient = column name.
raw_id  <- colnames(pb)
study_c <- as.character(pmeta[[opt$study_col]])
patient <- raw_id
is_ney  <- study_c == "neyazi"
patient[is_ney] <- vapply(strsplit(raw_id[is_ney], "_"), function(x) x[2], character(1))
n_rep <- sum(duplicated(patient))
message(sprintf("[DE] patient-level block: %d samples -> %d distinct patients (%d repeated-measure columns)",
                length(patient), length(unique(patient)), n_rep))

v <- voom(dge, design)
block <- patient
if (any(duplicated(block))) {
  corfit <- duplicateCorrelation(v, design, block = block)
  message(sprintf("[DE] duplicateCorrelation consensus = %.4f", corfit$consensus))
  fit <- lmFit(v, design, block = block, correlation = corfit$consensus)
} else {
  message("[DE] no repeated patients -> plain lmFit")
  fit <- lmFit(v, design)
}

lvl <- colnames(design)
mk  <- function(a, b) {
  ca <- paste0("disease", make.names(a)); cb <- paste0("disease", make.names(b))
  if (!(ca %in% lvl) || !(cb %in% lvl)) return(NULL)
  setNames((lvl == ca) - (lvl == cb), NULL)
}
contrasts <- list(
  CS_vs_DCM  = mk("CS", "DCM"),
  CS_vs_ARVC = mk("CS", "ARVC"),
  CS_vs_HCM  = mk("CS", "HCM"),
  CS_vs_NF   = mk("CS", "NF")
)
contrasts <- contrasts[!sapply(contrasts, is.null)]

for (nm in names(contrasts)) {
  cm_fit <- contrasts.fit(fit, contrasts[[nm]])
  cm_fit <- eBayes(cm_fit)
  tt <- topTable(cm_fit, number = Inf, sort.by = "P")
  contrast_confounded <- design_status != "full_rank" &&
    any(sapply(single_study, function(d) grepl(d, nm, ignore.case = TRUE)))
  tt$design_status <- design_status
  tt$contrast_confounded <- contrast_confounded
  out <- file.path(opt$out, sprintf("DE_%s.tsv", nm))
  write.table(tt, out, sep = "\t", quote = FALSE, col.names = NA)
  message(sprintf("wrote %s (%d genes, %d FDR<0.05)", out, nrow(tt),
                  sum(tt$adj.P.Val < 0.05, na.rm = TRUE)))
}
message("Patient-blocked DE complete.")
