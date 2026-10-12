#!/usr/bin/env Rscript
# Section-level pseudobulk limma-voom for the spatial-first analysis.
#   Q1: CS vs comparators, fitted on the 17 diseased sections.
# Detectability filter: detection fraction >= thresh in >= 1 of the two contrasted groups.
suppressPackageStartupMessages({library(limma); library(edgeR)})
a <- commandArgs(TRUE)
get <- function(k, d=NULL) { i <- which(a == k); if (length(i)) a[i+1] else d }
pb_f <- get("--pb"); det_f <- get("--det"); meta_f <- get("--meta")
thresh <- as.numeric(get("--thresh")); out <- get("--out"); tag <- get("--tag", "")
dir.create(out, recursive=TRUE, showWarnings=FALSE)

pb <- as.matrix(read.delim(pb_f, row.names=1, check.names=FALSE))   # sections x genes
meta <- read.delim(meta_f, row.names=1, check.names=FALSE)
det <- read.delim(det_f, row.names=1, check.names=FALSE)            # genes x groups

meta <- meta[meta$evaluable == "True" | meta$evaluable == TRUE, , drop=FALSE]
pb <- t(pb[rownames(meta), , drop=FALSE])                           # genes x sections
stopifnot(all(colnames(pb) == rownames(meta)))
cat(sprintf("[in] %d genes x %d evaluable sections\n", nrow(pb), ncol(pb)))
print(table(meta$cohort))

run_fit <- function(counts, m, design, contrasts, keep) {
  d <- DGEList(counts[keep, , drop=FALSE])
  d <- calcNormFactors(d, method="TMM")
  v <- voom(d, design, plot=FALSE)
  fit <- eBayes(contrasts.fit(lmFit(v, design), contrasts))
  list(fit=fit, n=sum(keep), v=v)
}

tt <- function(fit, coef) {
  t <- topTable(fit, coef=coef, number=Inf, sort.by="none")
  data.frame(gene=rownames(t), logFC=t$logFC, AveExpr=t$AveExpr, t=t$t,
             P.Value=t$P.Value, adj.P.Val=t$adj.P.Val, B=t$B)
}
# ---------------- Q1: CS vs comparators, diseased sections only ----------------
sel <- meta$cohort %in% c("CS", "COMPARATOR")
m1 <- meta[sel, , drop=FALSE]; pb1 <- pb[, sel, drop=FALSE]
coh <- factor(m1$cohort, levels=c("CS", "COMPARATOR"))
des1 <- model.matrix(~0 + coh); colnames(des1) <- levels(coh)
cont1 <- limma::makeContrasts(cs_vs_comp = CS - COMPARATOR, levels=des1)
keep1 <- (det$CS >= thresh) | (det$COMPARATOR >= thresh)
keep1[is.na(keep1)] <- FALSE
names(keep1) <- rownames(det); keep1 <- keep1[rownames(pb1)]
r1 <- run_fit(pb1, m1, des1, cont1, keep1)
cat(sprintf("[Q1] %d genes pass detectability >= %.2f\n", r1$n, thresh))
q1 <- tt(r1$fit, "cs_vs_comp")
ci1 <- 1.96 * sqrt(r1$fit$s2.post) * r1$fit$stdev.unscaled[, "cs_vs_comp"]
q1$CI_lo <- q1$logFC - ci1; q1$CI_hi <- q1$logFC + ci1
q1$sigma_post <- sqrt(r1$fit$s2.post)
q1$stdev_unscaled <- r1$fit$stdev.unscaled[, "cs_vs_comp"]
q1$df_total <- r1$fit$df.total
q1$det_CS <- det[q1$gene, "CS"]; q1$det_COMPARATOR <- det[q1$gene, "COMPARATOR"]
write.table(q1, file.path(out, paste0("q1_raw", tag, ".tsv")), sep="\t", quote=FALSE, row.names=FALSE)

cat(sprintf("[Q1] FDR<0.05: %d (up in CS %d)\n", sum(q1$adj.P.Val < 0.05),
            sum(q1$adj.P.Val < 0.05 & q1$logFC > 0)))
cat("[done]\n")
