#!/usr/bin/env Rscript
# ESC HF A5b: procurement in the snRNA arm. Run from repo root (env with limma + edgeR).
#
# Model (identical to the committed snRNA contrasts; reproduction asserted below):
#   patient-level CM pseudobulk counts (data/snrna/snrna_cm_patient_counts.tsv.gz),
#   filterByExpr(group = disease), TMM, voom, lmFit(~0 + disease + anatomy3), contrasts
#   CS - DCM, CS - ARVC, CS - NF, eBayes.
# (1) Transplant-only refit: CS restricted to Neyazi "Tissue Source" == TX; all comparator
#     patients kept; the same model re-run from filterByExpr onward. Pooling of DCM and ARVC
#     is done in a5b_pool.py with the 04_snrna_pool.py procedure.
# (2)/(3) per-patient scores are written here (log1p CP10K of the summed patient counts;
#     program = mean over QC-passing CS-lower genes present) and tested in a5b_pool.py.
suppressMessages({library(limma); library(edgeR)})
v2 <- Sys.getenv("ESCHF_PRIMARY", "old") == "v2"   # see eschf_paths.py
out <- if (v2) "results/eschf/primary_v2" else "results/eschf"
cnt <- as.matrix(read.delim(gzfile("data/snrna/snrna_cm_patient_counts.tsv.gz"), row.names = 1, check.names = FALSE))
coh <- read.delim("data/snrna/snrna_patient_cohort.tsv")
stopifnot(all(coh$patient %in% colnames(cnt)))
cnt <- cnt[, coh$patient]
clin <- read.delim("data/external/neyazi_clinical.tsv", check.names = FALSE)
coh$procurement <- clin$`Tissue Source`[match(coh$patient, clin$patient)]
stopifnot(all(!is.na(coh$procurement[coh$disease == "CS"])))

fit_contrasts <- function(keep_pat, label) {
  cc <- coh[coh$patient %in% keep_pat, ]
  y <- cnt[, cc$patient]
  disease <- factor(cc$disease); anatomy3 <- factor(cc$anatomy3)
  d <- DGEList(y)
  d <- d[filterByExpr(d, group = disease), , keep.lib.sizes = FALSE]
  d <- calcNormFactors(d)
  X <- model.matrix(~0 + disease + anatomy3)
  colnames(X) <- make.names(colnames(X))
  if (qr(X)$rank < ncol(X)) stop(sprintf("[%s] design rank %d < %d columns", label, qr(X)$rank, ncol(X)))
  v <- voom(d, X)
  f <- lmFit(v, X)
  cm <- makeContrasts(DCM = diseaseCS - diseaseDCM, ARVC = diseaseCS - diseaseARVC, NF = diseaseCS - diseaseNF,
                      levels = X)
  f2 <- eBayes(contrasts.fit(f, cm))
  res <- lapply(colnames(cm), function(k) {
    t <- topTable(f2, coef = k, number = Inf, sort.by = "none"); t$gene <- rownames(t); t })
  names(res) <- colnames(cm)
  cat(sprintf("[%s] patients=%d (CS=%d) genes=%d design=%dx%d rank=%d\n", label, ncol(y), sum(cc$disease == "CS"),
              nrow(d), nrow(X), ncol(X), qr(X)$rank))
  res
}

full <- fit_contrasts(coh$patient, "full")
for (k in c("DCM", "ARVC", "NF")) {
  ref <- read.delim(sprintf("data/snrna/snrna_CS_vs_%s.tsv", k))
  m <- match(ref$gene, full[[k]]$gene)
  stopifnot(!anyNA(m), max(abs(ref$logFC - full[[k]]$logFC[m])) < 1e-8, max(abs(ref$t - full[[k]]$t[m])) < 1e-8)
}
cat("[full] reproduces committed snrna_CS_vs_{DCM,ARVC,NF} (max |dlogFC|, |dt| < 1e-8)\n")

tx_keep <- coh$patient[coh$disease != "CS" | coh$procurement == "TX"]
tx <- fit_contrasts(tx_keep, "transplant_only")
for (k in c("DCM", "ARVC", "NF"))
  write.table(tx[[k]], file.path(out, sprintf("a5b_snrna_txonly_CS_vs_%s.tsv", k)), sep = "\t", quote = FALSE, row.names = FALSE)

# per-patient log1p CP10K scores (CS patients)
prog <- grep("^#", readLines(if (v2) file.path(out, "geneset_primary_program_qcpass_down.txt") else "metadata/geneset_primary_program_qcpass_down.txt"), invert = TRUE, value = TRUE)
prog <- prog[nzchar(prog)]
cs <- coh[coh$disease == "CS", ]
lib <- colSums(cnt[, cs$patient])
lc <- log1p(sweep(cnt[, cs$patient], 2, lib, "/") * 1e4)
pg <- intersect(prog, rownames(lc))
miss <- setdiff(prog, rownames(lc))
sc <- data.frame(patient = cs$patient, procurement = cs$procurement, n_cells = cs$n_cells,
                 program_score = colMeans(lc[pg, ]), ACADVL = lc["ACADVL", ], IDH2 = lc["IDH2", ])
sc$lvef_raw <- clin$lvef_raw[match(sc$patient, clin$patient)]
sc$lvef <- clin$lvef[match(sc$patient, clin$patient)]
write.table(sc, file.path(out, "a5b_cs_patient_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(analysis = "A5b snRNA program score", gene = miss, reason = "absent from snRNA count matrix"),
            file.path(out, "a5b_missing_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("[scores] program genes used %d of %d; procurement: %s\n", length(pg), length(prog),
            paste(names(table(sc$procurement)), table(sc$procurement), collapse = ", ")))
