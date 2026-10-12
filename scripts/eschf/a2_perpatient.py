#!/usr/bin/env python3
"""ESC HF revision s.3-4: per-patient distance slopes and the OXPHOS level check (new primary).

Run from the repo root after ESCHF_PRIMARY=v2 a2a6_analysis.py.
Specification (fixed in handoff_figure_fixes s.3, before results):
  For each CS patient with cardiomyocyte-dominant spots in >= 2 populated distance bins, OLS on that
  patient's spots:  score ~ dist_mm + cm_structural + z(log10 UMI) [+ section indicator if > 1 section].
  Positive control cm_structural: same model without the structural covariate.
  Outcomes: program score (new gene list), ACADVL, KEGG OXPHOS, CM-structural.
  Per outcome: n patients, mean per-patient slope with t-based 95% CI, counts of positive/negative slopes.
  Gap share (program, ACADVL): |mean slope x max distance| / (comparator mean - CS pooled mean), where the
  max distance is the largest CS cardiomyocyte-dominant spot distance in the data, the comparator mean is the
  mean of comparator patient means and the CS pooled mean is the mean of CS patient means over all
  cardiomyocyte-dominant spots (all bins pooled).
  z(log10 UMI) is standardised over all A2 spots, as in a2_model.R (a linear rescaling within a patient,
  so it does not change the distance slope).
OXPHOS level (s.4): per populated bin, CS mean of patient means with 95% t-interval beside the comparator
mean and range of patient means; whether the CS bin mean (and the upper CI) falls below the comparator minimum.
Outputs: results/eschf/a2_perpatient_slopes.tsv, a2_perpatient_slopes_bypatient.tsv, a2_oxphos_level.tsv
"""
import numpy as np
import pandas as pd
from scipy import stats

V2, OUT = "results/eschf/primary_v2", "results/eschf"
d = pd.read_csv(f"{V2}/a2_model_input.tsv.gz", sep="\t")
assert d.dist_gran_mm.notna().all(), "spots without a distance in the A2 input"
d["z_log10_umi"] = (d.log10_umi - d.log10_umi.mean()) / d.log10_umi.std(ddof=1)
d["cm_structural"] = d.score_cm_structural
OUTCOMES = {"program": ("score_program", True), "ACADVL": ("score_ACADVL", True),
            "oxphos_kegg": ("score_oxphos_kegg", True), "cm_structural_positive_control": ("score_cm_structural", False)}
bins_pop = d.groupby("patient").bin.nunique()
elig = sorted(bins_pop.index[bins_pop >= 2])

rows = []
for p in sorted(d.patient.unique()):
    sub = d[d.patient == p]
    for name, (col, with_struct) in OUTCOMES.items():
        r = dict(patient=p, outcome=name, n_spots=len(sub), n_sections=sub.section.nunique(),
                 n_populated_bins=int(bins_pop[p]), eligible=p in elig)
        if p in elig:
            X = [np.ones(len(sub)), sub.dist_gran_mm.values]
            if with_struct:
                X.append(sub.cm_structural.values)
            X.append(sub.z_log10_umi.values)
            for s in sorted(sub.section.unique())[1:]:
                X.append((sub.section == s).values.astype(float))
            X = np.column_stack(X)
            assert np.linalg.matrix_rank(X) == X.shape[1], f"rank-deficient design: {p} {name}"
            beta, *_ = np.linalg.lstsq(X, sub[col].values, rcond=None)
            r["slope_per_mm"] = beta[1]
            r["n_params"] = X.shape[1]
        rows.append(r)
pp = pd.DataFrame(rows)
pp.to_csv(f"{OUT}/a2_perpatient_slopes_bypatient.tsv", sep="\t", index=False)

cref = pd.read_csv(f"{V2}/a2_comparator_reference.tsv", sep="\t").set_index("score")
cs_pooled = d.groupby("patient")[[c for c, _ in OUTCOMES.values()]].mean().mean()
dmax = d.dist_gran_mm.max()
summ = []
for name, (col, _) in OUTCOMES.items():
    v = pp[(pp.outcome == name) & pp.eligible].slope_per_mm.values
    n = len(v); m = v.mean(); h = stats.t.ppf(0.975, n - 1) * v.std(ddof=1) / np.sqrt(n)
    r = dict(outcome=name, n_patients=n, mean_slope_per_mm=m, ci_low=m - h, ci_high=m + h,
             n_negative=int((v < 0).sum()), n_positive=int((v > 0).sum()),
             t_p=stats.ttest_1samp(v, 0).pvalue, max_distance_mm=dmax,
             comparator_mean=cref.loc[col, "mean_across_patients"], cs_pooled_mean=cs_pooled[col],
             n_cs_patients_pooled=d.patient.nunique())
    r["gap_comparator_minus_cs"] = r["comparator_mean"] - r["cs_pooled_mean"]
    r["gap_share_across_max_distance"] = (abs(m * dmax) / r["gap_comparator_minus_cs"]
                                          if name in ("program", "ACADVL") else np.nan)
    summ.append(r)
S = pd.DataFrame(summ)
S.to_csv(f"{OUT}/a2_perpatient_slopes.tsv", sep="\t", index=False)

b = pd.read_csv(f"{V2}/a2_distance_bins.tsv", sep="\t")
b = b[b.n_spots > 0]
c = cref.loc["score_oxphos_kegg"]
ox = pd.DataFrame(dict(bin=b.bin, n_spots=b.n_spots, n_patients=b.n_patients, n_sections=b.n_sections,
                       cs_mean_of_patient_means=b.score_oxphos_kegg_mean, cs_ci_low=b.score_oxphos_kegg_ci_low,
                       cs_ci_high=b.score_oxphos_kegg_ci_high, comparator_mean=c.mean_across_patients,
                       comparator_min=c["min"], comparator_max=c["max"], n_comparator_patients=int(c.n_patients)))
ox["cs_mean_below_comparator_min"] = ox.cs_mean_of_patient_means < ox.comparator_min
ox["cs_ci_high_below_comparator_min"] = ox.cs_ci_high < ox.comparator_min
ox.to_csv(f"{OUT}/a2_oxphos_level.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
print("eligible:", elig); print(S.round(4).T.to_string()); print(ox.round(4).to_string())
print(pp.pivot(index="patient", columns="outcome", values="slope_per_mm").round(4).to_string())
