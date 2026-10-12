#!/usr/bin/env Rscript
# ESC HF figure v3, A10 input: patient-pseudobulk log2 CPM (voom form) for the primary fit.
# Patient pseudobulk counts rebuilt exactly as in a0_primary_refit.R (retained sections summed within patient,
# tested universe, as in the patient refit of a0_primary_refit.R); TMM over the 14 fitted patients;
# log2((y + 0.5) / (N_eff + 1) * 1e6). Also asserts that cpm(prior.count = 2) reproduces a1_patient_logcpm.tsv.
suppressMessages(library(edgeR))
pb <- as.matrix(read.delim(gzfile("data/eschf/primary_v2/pseudobulk_allcm_T2000.tsv.gz"), row.names=1, check.names=FALSE))
q1 <- read.delim("results/eschf/primary_v2/q1_raw_allcm.tsv")
mp <- read.delim("results/eschf/primary_v2/section_patient_map.tsv")
ret <- mp[mp$retained == "yes", ]
Y0 <- t(pb[ret$section, q1$gene])
pts <- unique(ret$patient)
Y <- sapply(pts, function(p) rowSums(Y0[, ret$section[ret$patient == p], drop = FALSE]))
d <- calcNormFactors(DGEList(Y), method = "TMM")
ref <- as.matrix(read.delim("results/eschf/primary_v2/a1_patient_logcpm.tsv", row.names = 1, check.names = FALSE))
stopifnot(max(abs(cpm(d, log = TRUE, prior.count = 2) - ref[rownames(Y), colnames(Y)])) < 1e-8)
N <- d$samples$lib.size * d$samples$norm.factors
lc <- log2(t((t(Y) + 0.5) / (N + 1)) * 1e6)
write.table(data.frame(gene = rownames(lc), lc, check.names = FALSE), "results/eschf/a10_patient_log2cpm_voomform.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(dim(lc), "\n")
