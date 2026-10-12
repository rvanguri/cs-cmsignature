#!/usr/bin/env python3
"""ESC HF A2 (distance curve: bins, comparator reference) and A6 (tissue-level fibrosis).

Input: data/eschf/eschf_perspot_T2000.tsv.gz from a2a6_spot_export.py (cluster; verified to
re-derive the committed pseudobulk_T2000 and region counts exactly). Every spot there has
>= 2,000 native UMI and was thinned to 2,000 UMI; scores are mean log1p CP10K over member
genes on the thinned counts. Lineage call = argmax of the five pooled lineage z-scores.
Granuloma = immune-dominant spot (Myeloid or Lymphoid argmax) in a CS section
(user decision for A2). dist_gran_mm = distance to the nearest such spot in the same section.

A2 spots: every cardiomyocyte-dominant spot in every CS section, no immune-adjacency
exclusion. Bins fixed in advance. The spot-level mixed model is fitted in a2_model.R.
A6 spots: all spots, no cardiomyocyte restriction. Section -> patient by mean over sections.
CS vs comparator patients: difference in means with Welch 95% CI, Mann-Whitney U.
"""
import numpy as np
import pandas as pd
from scipy import stats

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from eschf_paths import OUT, V2  # noqa: E402
# v2: per-spot table with the new program / OXPHOS scores (a0_perspot_scores.py)
PERSPOT = "data/eschf/primary_v2/eschf_perspot_T2000.tsv.gz" if V2 else "data/eschf/eschf_perspot_T2000.tsv.gz"
BINS = [0, 0.25, 0.5, 1, 2, 3, np.inf]
LABELS = ["<0.25", "0.25-0.5", "0.5-1", "1-2", "2-3", ">3"]
SCORES = ["score_program", "score_ACADVL", "score_oxphos_kegg", "score_cm_structural", "score_IDH2"]

sp = pd.read_csv(PERSPOT, sep="\t")
pmap = pd.read_csv(f"{OUT}/section_patient_map.tsv", sep="\t").set_index("section").patient
sp["patient"] = sp.section.map(pmap)
assert sp.patient.notna().all()

# ------------------------------------------------------------------ A2
cs = sp[(sp.cohort == "CS") & sp.cm_dominant].copy()
has_gran = sp[sp.cohort == "CS"].groupby("section").immune_dominant.sum() > 0
cs["bin"] = pd.cut(cs.dist_gran_mm, BINS, right=False, labels=LABELS).astype(object)
cs.loc[~cs.section.map(has_gran), "bin"] = "no granuloma in section"
cs.loc[~cs.section.map(has_gran), "dist_gran_mm"] = np.nan
assert cs.bin.notna().all()
cs["log10_umi"] = np.log10(cs.umi_native)
cs[["spot", "section", "patient", "umi_native", "log10_umi", "dist_gran_mm", "bin"] + SCORES].to_csv(
    f"{OUT}/a2_model_input.tsv.gz", sep="\t", index=False)

order = LABELS + ["no granuloma in section"]
pp = cs.groupby(["bin", "patient"])[SCORES].mean()
pp["n_spots"] = cs.groupby(["bin", "patient"]).size()
pp = pp.reset_index()
pp.to_csv(f"{OUT}/a2_distance_bins_perpatient.tsv", sep="\t", index=False)
rows = []
for b in order:
    sub = pp[pp.bin == b]
    r = dict(bin=b, n_spots=int(sub.n_spots.sum()), n_patients=len(sub),
             n_sections=int(cs[cs.bin == b].section.nunique()))
    for s in SCORES:
        v = sub[s].values
        r[f"{s}_mean"] = v.mean() if len(v) else np.nan
        if len(v) >= 2:
            h = stats.t.ppf(0.975, len(v) - 1) * v.std(ddof=1) / np.sqrt(len(v))
            r[f"{s}_ci_low"], r[f"{s}_ci_high"] = v.mean() - h, v.mean() + h
        else:
            r[f"{s}_ci_low"] = r[f"{s}_ci_high"] = np.nan
        r[f"{s}_spot_mean"] = cs.loc[cs.bin == b, s].mean() if r["n_spots"] else np.nan
    rows.append(r)
bins = pd.DataFrame(rows)
bins.to_csv(f"{OUT}/a2_distance_bins.tsv", sep="\t", index=False)

sec = sp[sp.cohort == "CS"].groupby("section").agg(
    patient=("patient", "first"), n_spots=("spot", "size"), n_cm_dominant=("cm_dominant", "sum"),
    n_immune_dominant=("immune_dominant", "sum"), max_dist_gran_mm_all_spots=("dist_gran_mm", "max"))
sec["max_dist_gran_mm_cm_spots"] = cs.groupby("section").dist_gran_mm.max()
sec.to_csv(f"{OUT}/a2_section_max_distance.tsv", sep="\t")

cmp = sp[(sp.cohort == "COMPARATOR") & sp.cm_dominant]
cr = cmp.groupby("patient")[SCORES].mean()
cr["n_spots"] = cmp.groupby("patient").size()
cr["disease"] = cmp.groupby("patient").disease.first()
summ = pd.DataFrame({"mean_across_patients": cr[SCORES].mean(), "min": cr[SCORES].min(), "max": cr[SCORES].max()})
summ["n_patients"] = len(cr)
cr.reset_index().to_csv(f"{OUT}/a2_comparator_reference_perpatient.tsv", sep="\t", index=False)
summ.rename_axis("score").reset_index().to_csv(f"{OUT}/a2_comparator_reference.tsv", sep="\t", index=False)

# ------------------------------------------------------------------ A6
assert set(sp.cohort) == {"CS", "COMPARATOR"}, "normals are not in the export"
a = sp
s6 = a.groupby("section").agg(cohort=("cohort", "first"), patient=("patient", "first"), disease=("disease", "first"),
                              n_spots=("spot", "size"), frac_fibroblast_dominant=("fibroblast_dominant", "mean"),
                              mean_ecm7=("score_ecm7", "mean"), frac_immune_dominant=("immune_dominant", "mean"),
                              frac_cm_dominant=("cm_dominant", "mean"))
VARS = ["frac_fibroblast_dominant", "mean_ecm7", "frac_immune_dominant", "frac_cm_dominant"]
p6 = s6.groupby("patient").agg(cohort=("cohort", "first"), disease=("disease", "first"),
                               n_sections=("n_spots", "size"), n_spots=("n_spots", "sum"),
                               **{v: (v, "mean") for v in VARS})
res = []
for v in VARS:
    x = p6.loc[p6.cohort == "CS", v]; y = p6.loc[p6.cohort == "COMPARATOR", v]
    vx, vy = x.var(ddof=1) / len(x), y.var(ddof=1) / len(y)
    se = np.sqrt(vx + vy); df = (vx + vy) ** 2 / (vx ** 2 / (len(x) - 1) + vy ** 2 / (len(y) - 1))
    d = x.mean() - y.mean(); q = stats.t.ppf(0.975, df)
    res.append(dict(variable=v, n_CS_patients=len(x), n_comparator_patients=len(y), mean_CS=x.mean(),
                    median_CS=x.median(), mean_comparator=y.mean(), median_comparator=y.median(),
                    diff_CS_minus_comparator=d, ci_low=d - q * se, ci_high=d + q * se, welch_p=2 * stats.t.sf(abs(d / se), df),
                    mannwhitney_p=stats.mannwhitneyu(x, y, alternative="two-sided").pvalue))
R6 = pd.DataFrame(res)
with open(f"{OUT}/a6_tissue_fibrosis.tsv", "w") as fh:
    R6.to_csv(fh, sep="\t", index=False)
s6.to_csv(f"{OUT}/a6_tissue_fibrosis_persection.tsv", sep="\t")
p6.to_csv(f"{OUT}/a6_tissue_fibrosis_perpatient.tsv", sep="\t")
print(bins[["bin", "n_spots", "n_patients", "n_sections", "score_program_mean", "score_ACADVL_mean",
            "score_cm_structural_mean"]].to_string())
print(sec.to_string()); print(summ.to_string()); print(R6.to_string()); print(p6.round(3).to_string())
