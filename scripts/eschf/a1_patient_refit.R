#!/usr/bin/env Rscript
# ESC HF A1: sections -> patients, and patient-level refits of the primary CS-vs-comparator
# contrast. Run from the repo root.
#
# Gene universe is FIXED to the 985 genes tested in the original section-level fit
# (detectability >= 0.25 in CS or comparators, results/spatial_first/q1_raw_T2000_d25.tsv).
# Original fit is first reproduced from the committed pseudobulk and asserted identical.
#
# Refit 1, patient pseudobulk: counts summed over retained sections within patient,
#   TMM, voom, ~0 + arm.
# Refit 2, section level with patient blocking: TMM, voom -> duplicateCorrelation(block =
#   patient) -> voom(block, correlation) -> duplicateCorrelation again -> lmFit(block,
#   correlation)  (the two-pass procedure of the limma user's guide).
# QC-passing summaries restrict the 985-gene fit to QC-passing genes; FDR is not re-adjusted.
suppressPackageStartupMessages({library(limma); library(edgeR)})
out <- "results/eschf"
pb <- as.matrix(read.delim(gzfile("data/spatial_first/pseudobulk_T2000.tsv.gz"), row.names=1, check.names=FALSE))
meta <- read.delim("data/spatial_first/meta_T2000.tsv", check.names=FALSE)
rs <- read.delim("data/spatial_first/region_selection.tsv", check.names=FALSE)
rs <- rs[rs$T == 2000, ]
q1 <- read.delim("results/spatial_first/q1_raw_T2000_d25.tsv")
qc <- read.delim("results/spatial_first/q2_locus_qc.tsv")
hits <- read.delim(file.path(out, "q1_hits_annotated.tsv"))
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
m$granuloma_definition <- ifelse(m$cohort == "CS", "immune-dominant spots (lineage argmax)", "not applicable (no exclusion in comparators)")
m$has_granuloma <- ifelse(m$cohort == "CS", ifelse(m$n_granulomatous > 0, "yes", "no"), "NA")
m$n_immune_dominant_spots <- m$n_granulomatous
mp <- m[order(m$arm, m$section), c("section", "patient", "disease", "arm", "retained", "retained_spots",
                                    "has_granuloma", "n_immune_dominant_spots", "granuloma_definition")]
write.table(mp, file.path(out, "section_patient_map.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cnt <- do.call(rbind, lapply(c("CS", "COMPARATOR"), function(a) {
  x <- mp[mp$arm == a, ]
  data.frame(arm=a, sections_total=nrow(x), patients_total=length(unique(x$patient)),
             sections_retained=sum(x$retained == "yes"),
             patients_retained=length(unique(x$patient[x$retained == "yes"])),
             spots_retained=sum(x$retained_spots[x$retained == "yes"]))
}))
write.table(cnt, file.path(out, "a1_patient_counts.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
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
cat("[reproduce] original section-level fit reproduced exactly\n")

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
if (!is.finite(dc1$consensus.correlation)) {
  # refit 2 is not estimable (consensus correlation non-finite); recorded as NA, not substituted
  cat(sprintf("[dupcor] NOT ESTIMABLE: consensus correlation %s; patients with >1 section: %d\n",
              dc1$consensus.correlation, sum(table(blk) > 1)))
  dc2 <- list(consensus.correlation=NA_real_); f2 <- list(df.prior=NA_real_)
  o2 <- data.frame(gene=genes, logFC=NA_real_, t=NA_real_, P=NA_real_, adj=NA_real_, CI_lo=NA_real_, CI_hi=NA_real_)
} else {
v2 <- voom(d0, des0, block=blk, correlation=dc1$consensus.correlation)
dc2 <- duplicateCorrelation(v2, des0, block=blk)
f2 <- eBayes(contrasts.fit(lmFit(v2, des0, block=blk, correlation=dc2$consensus.correlation), C0))
o2 <- tidy(f2)
}
cat(sprintf("[dupcor] consensus correlation pass1 %.4f, pass2 %.4f\n", dc1$consensus.correlation, dc2$consensus.correlation))

# ---------------------------------------------------------------- summaries
qcpass <- genes %in% qc$gene[qc$qc_pass %in% c(TRUE, "True")]
names(qcpass) <- genes
hit192 <- q1$gene[q1$adj.P.Val < 0.05]
summ <- function(o, label, set) {
  o <- o[match(genes, o$gene), ]; ref <- q1[match(genes, q1$gene), ]
  s <- set
  h <- intersect(hit192, genes[s])
  data.frame(refit=label, gene_set=ifelse(all(s), "all_985", "qc_pass"),
             n_genes=sum(s), hits_fdr5=sum(o$adj[s] < 0.05), hits_up=sum(o$adj[s] < 0.05 & o$logFC[s] > 0),
             hits_down=sum(o$adj[s] < 0.05 & o$logFC[s] < 0),
             pearson_r_vs_original=cor(o$logFC[s], ref$logFC[s]),
             n_orig_hits=length(h),
             sign_concordant_orig_hits=sum(sign(o$logFC[match(h, o$gene)]) == sign(ref$logFC[match(h, ref$gene)])),
             n_orig_down_hits=sum(ref$logFC[match(h, ref$gene)] < 0),
             orig_down_hits_still_down=sum(o$logFC[match(h, o$gene)] < 0 & ref$logFC[match(h, ref$gene)] < 0),
             orig_hits_still_fdr5=sum(o$adj[match(h, o$gene)] < 0.05),
             orig_hits_ci_excludes_0=sum(o$CI_lo[match(h, o$gene)] > 0 | o$CI_hi[match(h, o$gene)] < 0))
}
all_s <- rep(TRUE, length(genes))
S <- rbind(summ(o0, "original_section_level", all_s), summ(o0, "original_section_level", qcpass),
           summ(o1, "patient_pseudobulk", all_s), summ(o1, "patient_pseudobulk", qcpass),
           summ(o2, "section_dupcor_patient_block", all_s), summ(o2, "section_dupcor_patient_block", qcpass))
S$n_units_CS <- c(rep(sum(arm == "CS"), 2), rep(sum(parm == "CS"), 2), rep(sum(arm == "CS"), 2))
S$n_units_COMPARATOR <- c(rep(sum(arm == "COMPARATOR"), 2), rep(sum(parm == "COMPARATOR"), 2), rep(sum(arm == "COMPARATOR"), 2))
S$dupcor_consensus <- c(NA, NA, NA, NA, dc2$consensus.correlation, dc2$consensus.correlation)
S$df_prior <- c(f0$df.prior, f0$df.prior, f1$df.prior, f1$df.prior, f2$df.prior, f2$df.prior)
write.table(S, file.path(out, "a1_patient_refit_summary.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
print(S)

pg <- data.frame(gene=genes, qc_pass=qcpass, original_hit=genes %in% hit192,
                 logFC_orig=q1$logFC, adj_orig=q1$adj.P.Val)
for (nm in c("patient", "dupcor")) {
  o <- if (nm == "patient") o1 else o2
  o <- o[match(genes, o$gene), ]
  pg[[paste0("logFC_", nm)]] <- o$logFC; pg[[paste0("CI_lo_", nm)]] <- o$CI_lo
  pg[[paste0("CI_hi_", nm)]] <- o$CI_hi; pg[[paste0("P_", nm)]] <- o$P; pg[[paste0("adj_", nm)]] <- o$adj
}
write.table(pg, file.path(out, "a1_patient_refit_pergene.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
# patient-level log-CPM for figure panel D (ACADVL per patient)
lc <- cpm(d1, log=TRUE, prior.count=2)
write.table(data.frame(gene=rownames(lc), lc, check.names=FALSE), file.path(out, "a1_patient_logcpm.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
cat("[done]\n")
