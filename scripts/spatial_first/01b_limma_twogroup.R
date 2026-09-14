#!/usr/bin/env Rscript
# Generic two-group section-level pseudobulk contrast (limma-voom, TMM), used for the
# internal technical controls: probe-panel-version contrast and chemistry-matched Q1 arms.
suppressPackageStartupMessages({library(limma); library(edgeR)})
a <- commandArgs(TRUE)
get <- function(k, d=NULL) { i <- which(a == k); if (length(i)) a[i+1] else d }
pb <- as.matrix(read.delim(get("--pb"), row.names=1, check.names=FALSE))
meta <- read.delim(get("--meta"), row.names=1, check.names=FALSE)
det <- read.delim(get("--det"), row.names=1, check.names=FALSE)
thresh <- as.numeric(get("--thresh")); g1 <- get("--g1"); g2 <- get("--g2")
detcols <- strsplit(get("--detcols"), ",")[[1]]
meta <- meta[!is.na(meta$group) & meta$group %in% c(g1, g2), , drop=FALSE]
pb <- t(pb[rownames(meta), , drop=FALSE])
grp <- factor(meta$group, levels=c(g1, g2))
covar <- get("--covar", NA)
if (!is.na(covar)) {
  z <- scale(as.numeric(meta[[covar]]))
  des <- model.matrix(~0 + grp + z); colnames(des)[1:2] <- levels(grp)
  cat(sprintf("[covar] adjusting for %s\n", covar))
} else {
  des <- model.matrix(~0 + grp); colnames(des) <- levels(grp)
}
cont <- limma::makeContrasts(contrasts = paste0(g1, "-", g2), levels = des)
keep <- rowSums(sapply(detcols, function(cc) det[rownames(pb), cc] >= thresh), na.rm=TRUE) > 0
d <- calcNormFactors(DGEList(pb[keep, , drop=FALSE]), method="TMM")
v <- voom(d, des, plot=FALSE)
fit <- eBayes(contrasts.fit(lmFit(v, des), cont))
t <- topTable(fit, number=Inf, sort.by="none")
out <- data.frame(gene=rownames(t), logFC=t$logFC, AveExpr=t$AveExpr, t=t$t,
                  P.Value=t$P.Value, adj.P.Val=t$adj.P.Val,
                  sigma_post=sqrt(fit$s2.post),
                  stdev_unscaled=fit$stdev.unscaled[,1], df_total=fit$df.total)
write.table(out, get("--outfile"), sep="\t", quote=FALSE, row.names=FALSE)
cat(sprintf("[%s vs %s] n1=%d n2=%d | %d genes tested | FDR<0.05: %d (up %d)\n",
            g1, g2, sum(grp==g1), sum(grp==g2), sum(keep),
            sum(out$adj.P.Val<0.05), sum(out$adj.P.Val<0.05 & out$logFC>0)))
