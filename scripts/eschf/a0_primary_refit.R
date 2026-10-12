#!/usr/bin/env Rscript
# ESC HF revision, section 0: NEW primary CS-vs-comparator contrast over ALL cardiomyocyte-dominant
# spots (no immune-adjacency exclusion; decision recorded before results, handoff_figure_fixes s.0).
# Run from the repo root.  Same code paths as a1_patient_refit.R; only the inputs change:
#   pseudobulk / detection / meta: data/eschf/primary_v2/*_allcm_T2000*  (a0_primary_export.py)
#   section-level fit:             results/eschf/primary_v2/q1_raw_allcm.tsv (01_limma_q1q2.R)
# Universe: detection >= 0.25 in CS or COMPARATOR CM-dominant spots (the committed rule).
# Fits: (0) section-level limma-voom (primary; reproduced here and asserted identical),
#       (1) patient pseudobulk, (2) section level with duplicateCorrelation(block = patient).
# QC: existing locus-QC flags (results/spatial_first/q2_locus_qc.tsv); genes outside that table are
#     not QC-passing (same convention as 00_inputs.py).  FDR is not re-adjusted within QC-passing.
suppressPackageStartupMessages({library(limma); library(edgeR)})
out <- "results/eschf"; V2DIR <- file.path(out, "primary_v2"); dir.create(V2DIR, showWarnings=FALSE)
pb <- as.matrix(read.delim(gzfile("data/eschf/primary_v2/pseudobulk_allcm_T2000.tsv.gz"), row.names=1, check.names=FALSE))
meta <- read.delim("data/eschf/primary_v2/meta_allcm_T2000.tsv", check.names=FALSE)
rs <- read.delim("data/spatial_first/region_selection.tsv", check.names=FALSE)
rs <- rs[rs$T == 2000, ]
q1 <- read.delim(file.path(V2DIR, "q1_raw_allcm.tsv"))
qc <- read.delim("results/spatial_first/q2_locus_qc.tsv")
genes <- q1$gene

# ---------------------------------------------------------------- section -> patient map
# CS: Foong sections named CS_<patient>[-<section>] (GSE314910: 12 sections, 9 patients).
# Comparators: six distinct patients (confirmed by the study PI, 2026-10-05).
m <- merge(meta[, c("section", "cohort", "disease", "n_selected", "evaluable")],
           rs[, c("section", "n_granulomatous")], by="section")
m <- m[m$cohort %in% c("CS", "COMPARATOR"), ]
m$patient <- ifelse(m$cohort == "CS", sub("-[0-9]+$", "", m$section), m$section)
m$arm <- m$cohort
m$retained <- ifelse(m$evaluable %in% c(TRUE, "True"), "yes", "no")
m$retained_spots <- m$n_selected
m$granuloma_definition <- ifelse(m$cohort == "CS", "immune-dominant spots (lineage argmax); not excluded in this fit", "not applicable")
m$has_granuloma <- ifelse(m$cohort == "CS", ifelse(m$n_granulomatous > 0, "yes", "no"), "NA")
m$n_immune_dominant_spots <- m$n_granulomatous
mp <- m[order(m$arm, m$section), c("section", "patient", "disease", "arm", "retained", "retained_spots",
                                    "has_granuloma", "n_immune_dominant_spots", "granuloma_definition")]
write.table(mp, file.path(V2DIR, "section_patient_map.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cnt <- do.call(rbind, lapply(c("CS", "COMPARATOR"), function(a) {
  x <- mp[mp$arm == a, ]
  data.frame(arm=a, sections_total=nrow(x), patients_total=length(unique(x$patient)),
             sections_retained=sum(x$retained == "yes"),
             patients_retained=length(unique(x$patient[x$retained == "yes"])),
             spots_retained=sum(x$retained_spots[x$retained == "yes"]))
}))
write.table(cnt, file.path(V2DIR, "a1_patient_counts.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
print(cnt)

ret <- mp[mp$retained == "yes", ]
arm <- factor(ret$arm, levels=c("CS", "COMPARATOR"))
Y <- t(pb[ret$section, genes])

check_rank <- function(des) {
  r <- qr(des)$rank
  if (r != ncol(des)) stop(sprintf("design rank %d < %d columns", r, ncol(des)))
  r
}
tidy <- function(fit, coef="d") {
  t <- topTable(fit, coef=coef, number=Inf, sort.by="none")
  data.frame(gene=rownames(t), logFC=t$logFC, t=t$t, P=t$P.Value, adj=t$adj.P.Val,
             CI_lo=t$logFC - qt(0.975, fit$df.total) * fit$stdev.unscaled[, coef] * sqrt(fit$s2.post),
             CI_hi=t$logFC + qt(0.975, fit$df.total) * fit$stdev.unscaled[, coef] * sqrt(fit$s2.post))
}

# ---------------------------------------------------------------- original (reproduction)
des0 <- model.matrix(~0 + arm); colnames(des0) <- levels(arm); check_rank(des0)
C0 <- makeContrasts(d = CS - COMPARATOR, levels=des0)
d0 <- calcNormFactors(DGEList(Y), method="TMM")
f0 <- eBayes(contrasts.fit(lmFit(voom(d0, des0), des0), C0))
o0 <- tidy(f0)
stopifnot(max(abs(o0$logFC - q1$logFC[match(o0$gene, q1$gene)])) < 1e-8)
cat("[reproduce] new section-level fit (01_limma_q1q2.R) reproduced exactly\n")

# ---------------------------------------------------------------- refit 1: patient pseudobulk
pts <- unique(ret$patient)
Yp <- sapply(pts, function(p) rowSums(Y[, ret$section[ret$patient == p], drop=FALSE]))
parm <- factor(ret$arm[match(pts, ret$patient)], levels=c("CS", "COMPARATOR"))
des1 <- model.matrix(~0 + parm); colnames(des1) <- levels(parm); check_rank(des1)
C1 <- makeContrasts(d = CS - COMPARATOR, levels=des1)
d1 <- calcNormFactors(DGEList(Yp), method="TMM")
f1 <- eBayes(contrasts.fit(lmFit(voom(d1, des1), des1), C1))
o1 <- tidy(f1)

# ---------------------------------------------------------------- refit 2: dupCor block = patient
blk <- ret$patient
v <- voom(d0, des0)
dc1 <- duplicateCorrelation(v, des0, block=blk)
v2 <- voom(d0, des0, block=blk, correlation=dc1$consensus.correlation)
dc2 <- duplicateCorrelation(v2, des0, block=blk)
f2 <- eBayes(contrasts.fit(lmFit(v2, des0, block=blk, correlation=dc2$consensus.correlation), C0))
o2 <- tidy(f2)
cat(sprintf("[dupcor] consensus correlation pass1 %.4f, pass2 %.4f\n", dc1$consensus.correlation, dc2$consensus.correlation))

# ---------------------------------------------------------------- summaries
qcpass <- genes %in% qc$gene[qc$qc_pass %in% c(TRUE, "True")]
names(qcpass) <- genes
hit_new <- q1$gene[q1$adj.P.Val < 0.05]
summ <- function(o, label, set) {
  o <- o[match(genes, o$gene), ]; ref <- q1[match(genes, q1$gene), ]
  s <- set
  h <- intersect(hit_new, genes[s])
  data.frame(refit=label, gene_set=ifelse(all(s), "all_tested", "qc_pass"),
             n_genes=sum(s), hits_fdr5=sum(o$adj[s] < 0.05), hits_up=sum(o$adj[s] < 0.05 & o$logFC[s] > 0),
             hits_down=sum(o$adj[s] < 0.05 & o$logFC[s] < 0),
             pearson_r_vs_section_level=cor(o$logFC[s], ref$logFC[s]),
             n_section_hits=length(h),
             sign_concordant_section_hits=sum(sign(o$logFC[match(h, o$gene)]) == sign(ref$logFC[match(h, ref$gene)])),
             n_section_down_hits=sum(ref$logFC[match(h, ref$gene)] < 0),
             section_down_hits_still_down=sum(o$logFC[match(h, o$gene)] < 0 & ref$logFC[match(h, ref$gene)] < 0),
             section_hits_still_fdr5=sum(o$adj[match(h, o$gene)] < 0.05),
             section_hits_ci_excludes_0=sum(o$CI_lo[match(h, o$gene)] > 0 | o$CI_hi[match(h, o$gene)] < 0))
}
all_s <- rep(TRUE, length(genes))
S <- rbind(summ(o0, "section_level_primary", all_s), summ(o0, "section_level_primary", qcpass),
           summ(o1, "patient_pseudobulk", all_s), summ(o1, "patient_pseudobulk", qcpass),
           summ(o2, "section_dupcor_patient_block", all_s), summ(o2, "section_dupcor_patient_block", qcpass))
S$n_units_CS <- c(rep(sum(arm == "CS"), 2), rep(sum(parm == "CS"), 2), rep(sum(arm == "CS"), 2))
S$n_units_COMPARATOR <- c(rep(sum(arm == "COMPARATOR"), 2), rep(sum(parm == "COMPARATOR"), 2), rep(sum(arm == "COMPARATOR"), 2))
S$dupcor_consensus <- c(NA, NA, NA, NA, dc2$consensus.correlation, dc2$consensus.correlation)
S$df_prior <- c(f0$df.prior, f0$df.prior, f1$df.prior, f1$df.prior, f2$df.prior, f2$df.prior)
write.table(S, file.path(out, "a0_primary_refit_summary.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
print(S)

pg <- data.frame(gene=genes, in_qc_universe=genes %in% qc$gene, qc_pass=qcpass,
                 logFC_section=q1$logFC, CI_lo_section=q1$CI_lo, CI_hi_section=q1$CI_hi,
                 P_section=q1$P.Value, adj_section=q1$adj.P.Val, AveExpr=q1$AveExpr,
                 det_CS=q1$det_CS, det_COMPARATOR=q1$det_COMPARATOR)
for (nm in c("patient", "dupcor")) {
  o <- if (nm == "patient") o1 else o2
  o <- o[match(genes, o$gene), ]
  pg[[paste0("logFC_", nm)]] <- o$logFC; pg[[paste0("CI_lo_", nm)]] <- o$CI_lo
  pg[[paste0("CI_hi_", nm)]] <- o$CI_hi; pg[[paste0("P_", nm)]] <- o$P; pg[[paste0("adj_", nm)]] <- o$adj
}
write.table(pg, file.path(out, "a0_primary_refit.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
# patient-level log-CPM for figure panel D (ACADVL per patient)
lc <- cpm(d1, log=TRUE, prior.count=2)
write.table(data.frame(gene=rownames(lc), lc, check.names=FALSE), file.path(V2DIR, "a1_patient_logcpm.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
cat("[done]\n")
