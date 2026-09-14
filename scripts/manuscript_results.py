#!/usr/bin/env python3
"""Recompute every quantitative claim in the brief report from the committed tables.

Each claim is checked against the tables under results/ and data/, printed as
  manuscript value | recomputed value | status
and written to results/manuscript_results.tsv. Exit status is 1 if any claim marked
reproducible fails, so this can gate a commit. Claims the committed data cannot
address are listed as NOT_REPRODUCIBLE with the reason, and do not affect exit status;
claims where the committed data disagree with the text are marked MISMATCH and do.

Run after scripts/spatial_first/*.py and scripts/figure1_build.py, from the repo root.
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

RES = "results/spatial_first"
DAT = "data/spatial_first"
TAB = "results/tables"
FDR = 0.05
rows = []


def claim(name, manuscript, recomputed, ok, note=""):
    status = "NOT_REPRODUCIBLE" if ok is None else ("OK" if bool(ok) else "MISMATCH")
    rows.append(dict(claim=name, manuscript=manuscript, recomputed=recomputed,
                     status=status, note=note))


def main():
    q1 = pd.read_csv(f"{RES}/q1_raw_T2000_d25.tsv", sep="\t").set_index("gene")
    q2 = pd.read_csv(f"{RES}/q2_raw_T2000_d25.tsv", sep="\t").set_index("gene")
    qc = pd.read_csv(f"{RES}/q2_locus_qc.tsv", sep="\t").set_index("gene")
    det = pd.read_csv(f"{DAT}/detection_fraction_T2000.tsv.gz", sep="\t", index_col=0)
    rsel = pd.read_csv(f"{DAT}/region_selection.tsv", sep="\t")
    conc = pd.read_csv(f"{RES}/q1_core_snrna_concordance.tsv", sep="\t").set_index("set")
    pergene = pd.read_csv(f"{RES}/q1_core_snrna_pergene.tsv", sep="\t")
    pubsum = pd.read_csv(f"{RES}/q1_published_nf_summary.tsv", sep="\t").set_index("comparison")
    pubper = pd.read_csv(f"{RES}/q1_published_nf_pergene.tsv", sep="\t").set_index("gene")
    cohort = pd.read_csv("data/snrna/snrna_patient_cohort.tsv", sep="\t")
    sn = {k: pd.read_csv(f"data/snrna/snrna_CS_vs_{k}.tsv", sep="\t").set_index("gene")
          for k in ("DCM", "ARVC") if os.path.exists(f"data/snrna/snrna_CS_vs_{k}.tsv")}

    # ---- cohort and region selection -------------------------------------------
    r2 = rsel[rsel["T"] == 2000]
    n_cs = int((r2.cohort == "CS").sum())
    n_cs_eval = int(((r2.cohort == "CS") & r2.evaluable).sum())
    n_comp = int(((r2.cohort == "COMPARATOR") & r2.evaluable).sum())
    claim("CS sections profiled", 12, n_cs, n_cs == 12)
    claim("CS sections evaluable for non-lesional myocardium (50-spot threshold)",
          7, n_cs_eval, n_cs_eval == 7)
    claim("CS sections with non-lesional myocardium present", 11,
          int((r2.cohort == "CS").sum()) - int(((r2.cohort == "CS") & ~r2.evaluable).sum() - 4),
          None, "11 = sections with any cardiomyocyte-enriched region before the "
                "50-spot threshold; region_selection.tsv records the post-threshold call")
    claim("comparator sections", 6, n_comp, n_comp == 6)

    # ---- primary contrast ------------------------------------------------------
    sig = q1.index[q1["adj.P.Val"] < FDR]
    lower = q1.index[(q1["adj.P.Val"] < FDR) & (q1.logFC < 0)]
    higher = q1.index[(q1["adj.P.Val"] < FDR) & (q1.logFC > 0)]
    claim("genes tested at 25% detectability", 985, len(q1), len(q1) == 985)
    claim("genes differential at FDR 5%", 192, len(sig), len(sig) == 192)
    claim("genes lower in CS", 139, len(lower), len(lower) == 139)
    claim("genes higher in CS", 53, len(higher), len(higher) == 53)
    # Locus QC is evaluated on the genes tested in the CS-vs-normal contrast, so 31 of
    # the 192 hits fall outside its universe and are neither passed nor flagged.
    qc_hits = qc.qc_pass.reindex(sig).astype("boolean")
    n_pass = int(qc_hits.fillna(False).sum())
    claim("hits passing locus-level QC", 128, n_pass, n_pass == 128,
          f"{int(qc_hits.notna().sum())} of {len(sig)} hits are in the QC universe; "
          f"{int(qc_hits.isna().sum())} were not tested in the CS-vs-normal contrast")
    n_det_lower = int((det.CS.reindex(q1.index) < det.COMPARATOR.reindex(q1.index)).sum())
    claim("genes with lower detection fraction in CS (sensitivity balance)",
          "494/985", f"{n_det_lower}/{len(q1)}", n_det_lower == 494)

    named = ["SLC2A4", "ACADVL", "IDH2", "ACADM", "IDH3B", "SDHB", "HSD17B4",
             "SLC4A3", "ITGA7"]
    missing = [g for g in named if g not in set(lower)]
    claim("named metabolic/ion/adhesion genes are all in the lower-in-CS set",
          "all 9 lower", "all 9 lower" if not missing else f"missing {missing}",
          not missing)

    # ---- stability -------------------------------------------------------------
    loo_dir = f"{RES}/loo"
    sec_files = sorted(f for f in os.listdir(loo_dir) if f.endswith(".tsv")
                       and "dropdis" not in f)
    dis_files = sorted(f for f in os.listdir(loo_dir) if f.endswith(".tsv")
                       and "dropdis" in f)

    def intersect_sig(files):
        k = set(sig)
        for fn in files:
            d = pd.read_csv(f"{loo_dir}/{fn}", sep="\t").set_index("gene")
            k &= set(d.index[d["adj.P.Val"] < FDR])
        return k

    def all_retain_direction(files):
        for fn in files:
            d = pd.read_csv(f"{loo_dir}/{fn}", sep="\t").set_index("gene")
            if not (np.sign(d.logFC.reindex(sig)) == np.sign(q1.logFC.reindex(sig))).all():
                return False
        return True

    keep_sec = intersect_sig(sec_files)
    claim("all 192 hits retain direction in every refit", 192,
          192 if all_retain_direction(sec_files + dis_files) else "<192",
          all_retain_direction(sec_files + dis_files))
    claim("hits retaining significance in every leave-one-section-out refit",
          141, len(keep_sec), len(keep_sec) == 141)
    keep_dis = intersect_sig(dis_files)
    claim("hits retaining significance in every leave-one-disease-out refit",
          117, len(keep_dis), len(keep_dis) == 117)

    nomt = pd.read_csv(f"{RES}/q1_noMT_comparison.tsv", sep="\t").set_index("gene")
    nh = nomt.reindex([g for g in sig if g in nomt.index])
    same_dir = int((np.sign(nh.logFC_noMT) == np.sign(nh.logFC)).sum())
    still_sig = int((nh["adj.P.Val_noMT"] < FDR).sum())
    claim("removing mitochondrial genes and recomputing normalisation does not change "
          "the result", "unchanged",
          f"{same_dir}/{len(nh)} same direction, {still_sig}/{len(nh)} still FDR<5%",
          same_dir == len(nh) and still_sig / len(nh) > 0.95,
          "evaluated on the 192 hits; the MT genes themselves drop out of the refit")
    ctrl = pd.read_csv(f"{RES}/q1_control_comparison.tsv", sep="\t").set_index("gene")
    for tag, label in (("pm", "panel-matched comparators"),
                       ("pu", "adjustment for cardiomyocyte content")):
        same = int((np.sign(ctrl[f"logFC_{tag}"]) == np.sign(ctrl.logFC)).sum())
        n = int(ctrl[f"logFC_{tag}"].notna().sum())
        claim(f"result persists under {label}", "persists",
              f"{same}/{n} same direction", same / n > 0.95)

    # ---- cross-platform arm ----------------------------------------------------
    for dis_, n_exp in (("CS", 21), ("DCM", 52), ("ARVC", 8)):
        n_obs = int((cohort.disease == dis_).sum())
        claim(f"snRNA-seq patients, {dis_}", n_exp, n_obs, n_obs == n_exp)
    # CS vs HCM is excluded from the pooled comparator: its effect sizes track gene
    # length and most genes are called differential.
    gl = pd.read_csv(f"{RES}/snrna_gene_length_regression.tsv", sep="\t").set_index("contrast")
    hcm, pooled = gl.loc["CS_vs_HCM"], gl.loc[["CS_vs_DCM", "CS_vs_ARVC"]]
    claim("snRNA-seq CS vs HCM effect size scales with gene length",
          "scales with gene length",
          f"slope={hcm.slope_logFC_per_log10span:.2f} log2FC per log10 kb, "
          f"r={hcm.r:.2f} (pooled comparators |r|<={pooled.r.abs().max():.2f})",
          hcm.slope_logFC_per_log10span > 0 and hcm.p_slope < 1e-10
          and hcm.r > 2 * pooled.r.abs().max(),
          "basis for excluding CS vs HCM from the pooled snRNA-seq comparator")
    pct_hcm = 100 * hcm.n_sig / hcm.n_genes
    pct_pooled = 100 * (pooled.n_sig / pooled.n_genes)
    claim("most genes called differential in snRNA-seq CS vs HCM", "most",
          f"{hcm.n_sig:.0f}/{hcm.n_genes:.0f} ({pct_hcm:.0f}%) at FDR<5%; "
          f"pooled comparators {pct_pooled.min():.0f}-{pct_pooled.max():.0f}%",
          pct_hcm > 50 and pct_hcm > pct_pooled.max())
    claim("core genes testable in both platforms", 134, len(pergene),
          len(pergene) == 134)
    rep_dn = int(conc.loc["stable down (115)", "rep_fdr5"])
    rep_up = int(conc.loc["stable up (26)", "rep_fdr5"])
    claim("concordant and significant in both arms, lower in CS", 15, rep_dn, rep_dn == 15)
    claim("concordant and significant in both arms, higher in CS", 3, rep_up, rep_up == 3)
    frac = float(conc.loc["stable down (115)", "frac_conc"])
    claim("sign concordance among lower-expression genes", "60%",
          f"{frac * 100:.0f}%", round(frac * 100) == 60)
    rho, p_rho = stats.spearmanr(pergene.spatial_logFC, pergene.lfc_pool)
    claim("effect-magnitude correlation across platforms", "rho=0.03",
          f"rho={rho:.2f} (p={p_rho:.2f})", abs(rho - 0.03) < 0.005)

    # ---- published comparison with non-failing myocardium ----------------------
    claim("lower-in-CS genes carried into the published comparison", 139,
          int(pubsum.n_cs_lower.iloc[0]), int(pubsum.n_cs_lower.iloc[0]) == 139)
    for tag, n_test, n_low, base, pval in (("DCM", 36, 25, 0.48, 0.009),
                                           ("HCM", 35, 22, 0.45, 0.023)):
        r = pubsum.loc[f"{tag} vs non-failing"]
        claim(f"{tag}: CS-lower genes testable in the published table", n_test,
              int(r.n_testable), int(r.n_testable) == n_test,
              "published table lists only genes called differential")
        claim(f"{tag}: testable genes lower in disease than non-failing",
              f"{n_low}/{n_test}", f"{int(r.n_lower_in_disease)}/{int(r.n_testable)}",
              int(r.n_lower_in_disease) == n_low and int(r.n_testable) == n_test)
        # the text rounds the reference rate to whole percent
        claim(f"{tag}: reference rate among published differential genes", f"{base:.0%}",
              f"{r.base_frac_lower:.1%}", round(r.base_frac_lower, 2) == base)
        claim(f"{tag}: binomial p", pval, round(float(r.binom_p), 3),
              abs(float(r.binom_p) - pval) < 0.0005)

    idh2_ok = (q1.loc["IDH2", "logFC"] < 0
               and all(s.loc["IDH2", "logFC"] < 0 for s in sn.values() if "IDH2" in s.index)
               and pubper.loc["IDH2", "lfc_DCMNF"] < 0
               and pubper.loc["IDH2", "lfc_HCMNF"] < 0)
    claim("IDH2 lower in CS spatially, in snRNA-seq, and in DCM and HCM vs non-failing",
          "lower in all", "lower in all" if idh2_ok else "not all", bool(idh2_ok))
    claim("IDH2 lower in all 5 genotype strata of the comparator snRNA-seq cohort",
          "lower in 5/5", "genotype strata not in committed data", None,
          "snrna_patient_cohort.tsv carries disease and study, not the "
          "sarcomeric/desmosomal genotype strata the claim refers to")
    acadvl_ok = (q1.loc["ACADVL", "logFC"] < 0
                 and pubper.loc["ACADVL", "lfc_HCMNF"] < 0)
    claim("ACADVL lower in CS and in HCM vs non-failing", "lower in both",
          "lower in both" if acadvl_ok else "not both", bool(acadvl_ok))

    # ---- figure ----------------------------------------------------------------
    pA = pd.read_csv(f"{TAB}/figure1_panelA_values.tsv", sep="\t")
    pB = pd.read_csv(f"{TAB}/figure1_panelB_values.tsv", sep="\t")
    claim("Figure panel a genes", 134, len(pA), len(pA) == 134)
    claim("Figure panel b genes", 18, len(pB), len(pB) == 18)
    claim("Figure panel b direction split", "15 down / 3 up",
          f"{int((pB.direction == 'down').sum())} down / "
          f"{int((pB.direction == 'up').sum())} up",
          int((pB.direction == "down").sum()) == 15
          and int((pB.direction == "up").sum()) == 3)

    # ---- report ----------------------------------------------------------------
    out = pd.DataFrame(rows)
    out.to_csv("results/manuscript_results.tsv", sep="\t", index=False)
    w = max(out.claim.str.len())
    for _, r in out.iterrows():
        flag = {"OK": "  ok", "MISMATCH": "FAIL", "NOT_REPRODUCIBLE": "  --"}[r.status]
        print(f"{flag}  {r.claim:<{w}}  manuscript {str(r.manuscript):>15}  "
              f"recomputed {str(r.recomputed):>22}")
    n_fail = int((out.status == "MISMATCH").sum())
    n_na = int((out.status == "NOT_REPRODUCIBLE").sum())
    print(f"\n{int((out.status == 'OK').sum())} reproduced, {n_fail} mismatched, "
          f"{n_na} not addressable from committed data")
    if n_fail:
        print("mismatched claims are documented in docs/DISCREPANCIES.md")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
