#!/usr/bin/env Rscript
# ESC HF A2 model. Input: results/eschf/a2_model_input.tsv.gz (a2a6_analysis.py).
# CS cardiomyocyte-dominant spots, all CS sections with a granuloma (every section has one).
#   score ~ distance_mm + cm_structural + z(log10 UMI) + (1 | patient/section)
#   positive control: cm_structural ~ distance_mm + z(log10 UMI) + (1 | patient/section)
# REML, Satterthwaite df (lmerTest); 95% CI = beta +/- t(df) * SE.
suppressMessages({library(lme4); library(lmerTest)})
out <- if (Sys.getenv("ESCHF_PRIMARY", "old") == "v2") "results/eschf/primary_v2" else "results/eschf"  # eschf_paths.py
d <- read.delim(gzfile(file.path(out, "a2_model_input.tsv.gz")))
d <- d[!is.na(d$dist_gran_mm), ]
d$z_log10_umi <- as.numeric(scale(d$log10_umi))
d$cm_structural <- d$score_cm_structural
fits <- list(
  program        = score_program ~ dist_gran_mm + cm_structural + z_log10_umi + (1 | patient/section),
  ACADVL         = score_ACADVL ~ dist_gran_mm + cm_structural + z_log10_umi + (1 | patient/section),
  oxphos_kegg    = score_oxphos_kegg ~ dist_gran_mm + cm_structural + z_log10_umi + (1 | patient/section),
  cm_structural_positive_control = cm_structural ~ dist_gran_mm + z_log10_umi + (1 | patient/section))
rows <- lapply(names(fits), function(k) {
  m <- lmer(fits[[k]], data = d, REML = TRUE)
  X <- getME(m, "X"); stopifnot(qr(X)$rank == ncol(X))
  cf <- summary(m)$coefficients
  b <- cf["dist_gran_mm", ]
  q <- qt(0.975, b["df"])
  vc <- as.data.frame(VarCorr(m))
  data.frame(outcome = k, term = "dist_gran_mm (per mm)", beta = b["Estimate"], se = b["Std. Error"],
             df_satterthwaite = b["df"], ci_low = b["Estimate"] - q * b["Std. Error"],
             ci_high = b["Estimate"] + q * b["Std. Error"], p = b["Pr(>|t|)"],
             n_spots = nrow(d), n_sections = length(unique(d$section)), n_patients = length(unique(d$patient)),
             sd_outcome = sd(model.response(model.frame(m))),
             var_patient = vc$vcov[vc$grp == "patient"], var_section = vc$vcov[vc$grp == "section:patient"],
             var_resid = vc$vcov[vc$grp == "Residual"], singular = isSingular(m), row.names = NULL)
})
res <- do.call(rbind, rows)
res$range_mm <- paste(round(range(d$dist_gran_mm), 3), collapse = "-")
write.table(res, file.path(out, "a2_distance_model.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(res, digits = 4)
