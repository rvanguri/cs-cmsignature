#!/usr/bin/env Rscript
# Leave-one-out refits of the Q1 contrast (CS vs comparator cardiomyopathies).
#
# Two families of refit, both re-fitting the SAME model as 01_limma_q1q2.R on a reduced
# set of comparator sections:
#   --mode section  drop one comparator SECTION at a time   -> q1_raw_drop_<section>.tsv
#   --mode disease  drop one comparator DISEASE at a time   -> q1_raw_dropdis_<disease>.tsv
#                   (diseases contributing a single section are covered by --mode section,
#                    so only multi-section diseases get their own file)
#
# The detectability filter is computed ONCE on the full cohort and held fixed, so every
# refit tests the identical 985-gene set and gene sets are directly comparable.
suppressPackageStartupMessages({library(limma); library(edgeR)})
a <- commandArgs(TRUE)
get <- function(k, d=NULL) { i <- which(a == k); if (length(i)) a[i+1] else d }
pb_f <- get("--pb"); det_f <- get("--det"); meta_f <- get("--meta")
thresh <- as.numeric(get("--thresh", "0.25")); out <- get("--out"); mode <- get("--mode", "section")
dir.create(out, recursive=TRUE, showWarnings=FALSE)

pb <- as.matrix(read.delim(pb_f, row.names=1, check.names=FALSE))
meta <- read.delim(meta_f, row.names=1, check.names=FALSE)
det <- read.delim(det_f, row.names=1, check.names=FALSE)
meta <- meta[meta$evaluable == "True" | meta$evaluable == TRUE, , drop=FALSE]
pb <- t(pb[rownames(meta), , drop=FALSE])

keep <- (det$CS >= thresh) | (det$COMPARATOR >= thresh)
keep[is.na(keep)] <- FALSE
names(keep) <- rownames(det)

fit_q1 <- function(m1, pb1) {
  coh <- factor(m1$cohort, levels=c("CS", "COMPARATOR"))
  des <- model.matrix(~0 + coh); colnames(des) <- levels(coh)
  cont <- limma::makeContrasts(cs_vs_comp = CS - COMPARATOR, levels=des)
  k <- keep[rownames(pb1)]
  d <- DGEList(pb1[k, , drop=FALSE])
  d <- calcNormFactors(d, method="TMM")
  v <- voom(d, des, plot=FALSE)
  fit <- eBayes(contrasts.fit(lmFit(v, des), cont))
  t <- topTable(fit, coef="cs_vs_comp", number=Inf, sort.by="none")
  q <- data.frame(gene=rownames(t), logFC=t$logFC, AveExpr=t$AveExpr, t=t$t,
                  P.Value=t$P.Value, adj.P.Val=t$adj.P.Val, B=t$B)
  ci <- 1.96 * sqrt(fit$s2.post) * fit$stdev.unscaled[, "cs_vs_comp"]
  q$CI_lo <- q$logFC - ci; q$CI_hi <- q$logFC + ci
  q$sigma_post <- sqrt(fit$s2.post)
  q$stdev_unscaled <- fit$stdev.unscaled[, "cs_vs_comp"]
  q$df_total <- fit$df.total
  q$det_CS <- det[q$gene, "CS"]; q$det_COMPARATOR <- det[q$gene, "COMPARATOR"]
  q
}

dis1 <- meta[meta$cohort %in% c("CS", "COMPARATOR"), , drop=FALSE]
comp_sections <- rownames(dis1)[dis1$cohort == "COMPARATOR"]

if (mode == "section") {
  targets <- comp_sections
  drop_of <- function(x) x
  prefix <- "q1_raw_drop_"
} else {
  tab <- table(dis1$disease[dis1$cohort == "COMPARATOR"])
  targets <- names(tab)[tab > 1]   # single-section diseases are already covered by --mode section
  drop_of <- function(x) rownames(dis1)[dis1$cohort == "COMPARATOR" & dis1$disease == x]
  prefix <- "q1_raw_dropdis_"
}

for (tg in targets) {
  drop <- drop_of(tg)
  m1 <- dis1[!(rownames(dis1) %in% drop), , drop=FALSE]
  q <- fit_q1(m1, pb[, rownames(m1), drop=FALSE])
  f <- file.path(out, paste0(prefix, tg, ".tsv"))
  write.table(q, f, sep="\t", quote=FALSE, row.names=FALSE)
  cat(sprintf("[%s] dropped %-18s -> %d sections left, %d genes, FDR<0.05: %d (up in CS %d)\n",
              mode, tg, ncol(pb[, rownames(m1), drop=FALSE]), nrow(q),
              sum(q$adj.P.Val < 0.05), sum(q$adj.P.Val < 0.05 & q$logFC > 0)))
}
cat("[done]\n")
