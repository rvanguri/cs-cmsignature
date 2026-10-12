#!/usr/bin/env Rscript
# ESC HF figure v3, s.3b input: per-study mean log2 CPM of the snRNA cardiomyocyte pseudobulk.
# Pre-specified: patient-level CM pseudobulk counts (data/snrna/snrna_cm_patient_counts.tsv.gz, 121 patients),
# TMM normalisation over all 121 patients (edgeR calcNormFactors), log2 CPM in the voom form
# log2((y + 0.5) / (N_eff + 1) * 1e6) with N_eff = library size x TMM factor, no gene filtering.
# Per-study mean = mean over all patients of that study (all diseases). Also writes the summed
# CALM1+CALM2+CALM3 row (counts summed before log2 CPM).
suppressMessages(library(edgeR))
cnt <- as.matrix(read.delim(gzfile("data/snrna/snrna_cm_patient_counts.tsv.gz"), row.names = 1, check.names = FALSE))
coh <- read.delim("data/snrna/snrna_patient_cohort.tsv")
stopifnot(identical(colnames(cnt), coh$patient))
cnt <- rbind(cnt, CALM1_CALM2_CALM3_sum = colSums(cnt[c("CALM1", "CALM2", "CALM3"), ]))
dge <- calcNormFactors(DGEList(cnt[rownames(cnt) != "CALM1_CALM2_CALM3_sum", ]))
N <- dge$samples$lib.size * dge$samples$norm.factors
lc <- log2(t((t(cnt) + 0.5) / (N + 1)) * 1e6)
studies <- sort(unique(coh$study))
m <- sapply(studies, function(s) rowMeans(lc[, coh$study == s, drop = FALSE]))
n <- table(coh$study)[studies]
out <- data.frame(gene = rownames(m), m, check.names = FALSE)
write.table(out, gzfile("results/eschf/a9_study_mean_log2cpm.tsv.gz"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("genes", nrow(m), "studies", paste(studies, n, sep = ":"), "\n")
print(round(m[c("FPGT", "NDUFC2", "CALM2", "HSPB2", "CRYAB", "CALM1_CALM2_CALM3_sum"), ], 2))
