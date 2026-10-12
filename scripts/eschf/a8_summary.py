#!/usr/bin/env python3
"""Write results/eschf/SUMMARY.md from the ESC HF tables (no hand-typed numbers).

Order (HANDOFF_figure_fixes s.0): the revision items and the new primary contrast first (all
cardiomyocyte-dominant spots, tables in results/eschf/primary_v2 plus the a0_/a2_perpatient/a2_oxphos tables
at the top level), then the earlier >=0.5 mm fit as a sensitivity analysis (tables at the top level of
results/eschf). Run from the repo root after the A0-A7 scripts.
"""
import numpy as np
import pandas as pd

TOP, V2 = "results/eschf", "results/eschf/primary_v2"
FDR = 0.05
rdp = lambda d, f: pd.read_csv(f"{d}/{f}", sep="\t")
f2 = lambda x: f"{x:.2f}"
f3 = lambda x: f"{x:.3f}"
pc = lambda x: f"{100 * x:.0f}%"
truth = lambda s: s.astype(str).str.upper().eq("TRUE")


def pfmt(p):
    return "p<0.001" if p < 0.001 else f"p={p:.3f}"


def refit_summary(path):
    """Harmonise the two refit-summary layouts (old: *_orig_*; new: *_section_*)."""
    s = pd.read_csv(path, sep="\t")
    s = s.rename(columns={"n_section_hits": "n_orig_hits", "sign_concordant_section_hits": "sign_concordant_orig_hits",
                          "section_hits_still_fdr5": "orig_hits_still_fdr5",
                          "pearson_r_vs_section_level": "pearson_r_vs_original"})
    s["refit"] = s.refit.replace({"section_level_primary": "section_level", "original_section_level": "section_level"})
    s["gene_set"] = s.gene_set.replace({"all_985": "all_tested"})
    for col in ["n_genes", "hits_fdr5", "hits_up", "hits_down", "n_orig_hits", "sign_concordant_orig_hits",
                "orig_hits_still_fdr5", "n_units_CS", "n_units_COMPARATOR"]:
        s[col] = s[col].astype("Int64")  # a not-estimable refit is NA
    return s.set_index(["refit", "gene_set"])


def fit_sections(D, refit_path, h, label, n0, prog_path):
    """Sections A1, A2 (spot-level), A5b, A3/A4, A5 for the fit whose tables are in D. h: section-level hit table."""
    rd = lambda f: rdp(D, f)
    L = []
    hits = h[h["adj.P.Val"] < FDR]; hq = hits[truth(hits.qc_pass)]
    nconc = len(rd("a4_snrna_nf_concordant18.tsv"))
    # A1
    c = rd("a1_patient_counts.tsv").set_index("arm")
    s = refit_summary(refit_path)
    L += [f"### {n0}.1 A1: patient-level analysis ({label})", ""]
    L += [f"- CS: {c.loc['CS','sections_total']} sections from {c.loc['CS','patients_total']} patients; "
          f"{c.loc['CS','sections_retained']} sections from {c.loc['CS','patients_retained']} patients pass the "
          f"evaluability threshold ({c.loc['CS','spots_retained']} spots).",
          f"- Comparators: {c.loc['COMPARATOR','sections_retained']} sections from {c.loc['COMPARATOR','patients_retained']} "
          f"distinct patients ({c.loc['COMPARATOR','spots_retained']} spots).", ""]
    L += ["| Fit | Gene set | Units CS/comp | Hits FDR<5% | Section-level hits same sign | Section-level hits still FDR<5% | Pearson r vs section level |",
          "|---|---|---|---|---|---|---|"]
    for r in s.itertuples():
        rf, gs = r.Index
        L.append(f"| {rf} | {gs} | {r.n_units_CS}/{r.n_units_COMPARATOR} | {r.hits_fdr5} ({r.hits_down} down, {r.hits_up} up) | "
                 f"{r.sign_concordant_orig_hits}/{r.n_orig_hits} | {r.orig_hits_still_fdr5}/{r.n_orig_hits} | {r.pearson_r_vs_original:.3f} |")
    pp = next(s.loc[[("patient_pseudobulk", "all_tested")]].itertuples())
    dc = next(s.loc[[("section_dupcor_patient_block", "all_tested")]].itertuples())
    L += ["", f"At patient level (sections summed per patient), {pp.hits_fdr5} of {pp.n_genes} genes are FDR<5% against "
          f"{s.loc[('section_level','all_tested'),'hits_fdr5']} at section level. "
          f"{pp.orig_hits_still_fdr5} of the {pp.n_orig_hits} section-level hits remain FDR<5%, and "
          f"{pp.sign_concordant_orig_hits} keep their sign. With sections retained and patient blocking "
          f"(duplicateCorrelation, consensus ρ = {dc.dupcor_consensus:.3f}), {dc.hits_fdr5} are FDR<5%.", ""]
    # A2 spot-level
    b = rd("a2_distance_bins.tsv"); m = rd("a2_distance_model.tsv").set_index("outcome")
    cref = rd("a2_comparator_reference.tsv").set_index("score")
    pop = b[b.n_spots > 0]
    L += [f"### {n0}.2 A2: score by distance to the nearest immune-dominant spot ({label} program genes)", ""]
    pg = [l.strip() for l in open(prog_path) if l.strip() and not l.startswith("#")]
    L += [f"Program score ({len(pg)} genes): mean of patient means "
          f"{f2(pop.score_program_mean.min())}–{f2(pop.score_program_mean.max())} across the populated bins, against a comparator "
          f"range of {f2(cref.loc['score_program','min'])}–{f2(cref.loc['score_program','max'])}. *ACADVL*: "
          f"{f2(pop.score_ACADVL_mean.min())}–{f2(pop.score_ACADVL_mean.max())}, against "
          f"{f2(cref.loc['score_ACADVL','min'])}–{f2(cref.loc['score_ACADVL','max'])}.", ""]
    L += ["| Bin | Spots | Patients | Program score (95% CI) | ACADVL (95% CI) |", "|---|---|---|---|---|"]
    for _, r in pop.iterrows():
        L.append(f"| {r.bin} | {r.n_spots} | {r.n_patients} | {f3(r.score_program_mean)} ({f3(r.score_program_ci_low)}–{f3(r.score_program_ci_high)}) | "
                 f"{f3(r.score_ACADVL_mean)} ({f3(r.score_ACADVL_ci_low)}–{f3(r.score_ACADVL_ci_high)}) |")
    L += ["", "Spot-level model (singular fits; intervals anti-conservative): outcome ~ distance + CM-structural score + "
          "z(log10 UMI) + (1 | patient/section), REML, Satterthwaite df.", "",
          "| Outcome | β per mm | 95% CI | p | Satterthwaite df | Patient variance | Singular |", "|---|---|---|---|---|---|---|"]
    for k, r in m.iterrows():
        L.append(f"| {k} | {r.beta:+.4f} | {r.ci_low:+.4f} to {r.ci_high:+.4f} | {r.p:.2g} | {r.df_satterthwaite:.0f} | "
                 f"{r.var_patient:.2g} | {r.singular} |")
    L += ["", f"{int(truth(m.singular).sum())} of {len(m)} fits are singular. The model has no random slope, so the "
          "Satterthwaite df are spot-level.", ""]
    # A5b
    t = rd("a5b_txonly_concordance.tsv").set_index(["gene_set", "cohort"])
    pr = rd("a5b_procurement.tsv")
    L += [f"### {n0}.3 A5b: procurement (VAD vs transplant) in the snRNA CS arm ({label})", ""]
    L += ["| Gene set | Cohort | Testable | Same sign as spatial | Null p | Replicate FDR<5% | Same sign as full | Spearman full vs TX-only |",
          "|---|---|---|---|---|---|---|---|"]
    for (gs, co), r in t.iterrows():
        L.append(f"| {gs.replace('concordant_core', f'concordant_{nconc}')} | {co} | {r.n_testable} | {r.n_concordant_with_spatial} ({pc(r.frac_concordant)}) | "
                 f"{'' if pd.isna(r.null_p) else f'{r.null_p:.3f}'} | {r.n_replicate_fdr5} | "
                 f"{'' if pd.isna(r.n_same_sign_as_full) else int(r.n_same_sign_as_full)} | "
                 f"{'' if pd.isna(r.spearman_full_vs_tx) else f3(r.spearman_full_vs_tx)} |")
    L.append("")
    for _, r in pr[pr.analysis == "A5b.2 LVAD vs TX within CS"].iterrows():
        L.append(f"- LVAD − TX within CS, {r.variable}: {r.diff_LVAD_minus_TX:+.3f} (95% CI {r.ci_low:+.3f} to {r.ci_high:+.3f}; "
                 f"n = {int(r.n_LVAD)} LVAD, {int(r.n_TX)} TX; Welch {pfmt(r.p)}).")
    for _, r in pr[pr.analysis.str.startswith("A5b.3 LVEF")].iterrows():
        L.append(f"- LVEF vs {r.variable} across CS patients: Spearman ρ = {r.spearman_rho:+.3f} ({pfmt(r.p)}, n = {int(r.n)}).")
    c18t = t.loc[("concordant_core", "transplant_only")]
    L += ["- Disease duration: not available in the Neyazi supplementary tables.",
          f"- Transplant-only refit, {nconc} concordant genes: {int(c18t.n_same_sign_as_full)} of {nconc} keep their sign; "
          f"{int(c18t.n_replicate_fdr5)} of {nconc} remain FDR<5% "
          f"(full cohort: {int(t.loc[('concordant_core','full'),'n_replicate_fdr5'])}).", ""]
    # A3/A4
    a3 = rd("a3_qc_comparison.tsv")
    L += [f"### {n0}.4 A3 QC comparison and A4 non-failing ({label})", "",
          f"Hits are spatial FDR<5% ({int((hits.logFC < 0).sum())} lower, {int((hits.logFC > 0).sum())} higher in CS; "
          f"{int((hq.logFC < 0).sum())} and {int((hq.logFC > 0).sum())} pass locus QC). "
          "snRNA = inverse-variance pooled CS vs DCM and CS vs ARVC.", ""]
    keep = ["n_testable", "n_concordant", "frac_concordant", "null_frac", "null_p", "spearman_rho", "n_lower_in_disease",
            "frac_lower_in_disease", "ref_frac_lower_in_disease", "binom_p", "n_lower_in_CS_vs_NF", "frac_lower_in_CS_vs_NF",
            "ref_frac_lower_jointly_testable"]
    sub = a3[a3.statistic.isin(keep) & ~a3.analysis.str.startswith("hamed")]
    L += ["| Analysis | Gene universe | Direction | Statistic | All hits | QC-passing only |", "|---|---|---|---|---|---|"]
    for _, r in sub.iterrows():
        fmt = (lambda v: f"{int(v)}") if r.statistic.startswith("n_") else (lambda v: f"{v:.3g}")
        L.append(f"| {r.analysis} | {r.gene_universe} | {r.direction} | {r.statistic} | {fmt(r.all_hits)} | {fmt(r.qc_passing_only)} |")
    nf = rd("a4_internal_cs_vs_nf.tsv").set_index("statistic")
    L += ["", f"Internal snRNA CS vs NF: {int(nf.loc['n_lower_in_CS_vs_NF','qc_passing_only'])} of "
          f"{int(nf.loc['n_testable','qc_passing_only'])} QC-passing CS-lower genes are lower in CS vs NF "
          f"({pc(nf.loc['frac_lower_in_CS_vs_NF','qc_passing_only'])}; background {pc(nf.loc['ref_frac_lower_jointly_testable','qc_passing_only'])}).",
          f"Per-gene snRNA NF results for the {nconc} concordant genes are in `a4_snrna_nf_concordant18.tsv`; genes excluded "
          "from the QC-pass subset, with reasons, in `a3_excluded_genes.tsv`.", ""]
    # A5
    hv = rd("a5_hamed_vs_control.tsv"); sev = rd("a5_severity_axis.tsv")
    L += [f"### {n0}.5 A5 Hamed proteomics ({label})", "",
          "CS-lower spatial hits quantified in the proteome vs background (other tested genes in the proteome), "
          "fraction lower in each etiology vs controls (hypergeometric, both tails):", "",
          "| Etiology | QC | In proteome | Lower in disease | Background lower | p (depletion) | p (enrichment) |", "|---|---|---|---|---|---|---|"]
    for _, r in hv.iterrows():
        L.append(f"| {r.etiology} | {r.qc} | {r.n_in_proteome} | {r.n_lower_in_disease} ({pc(r.frac_lower)}) | {pc(r.bg_frac_lower)} | "
                 f"{r.p_depletion:.3f} | {r.p_enrichment:.3f} |")
    L += ["", "Severity axis (NMF Group 2 − Group 1, OLS with batch on unimputed values; CM-adjusted adds the sarcomeric content score):", "",
          "| Model | Test | QC | n | Value | Background | p (depletion) | p (enrichment) / p |", "|---|---|---|---|---|---|---|---|"]
    for _, r in sev.iterrows():
        pv = r.p if pd.notna(r.p) else r.p_enrichment
        L.append(f"| {r.model} | {r.test} | {r.qc} | {int(r.n)} | {r.value:.3g}{'' if pd.isna(r.frac) else f' ({pc(r.frac)})'} | "
                 f"{'' if pd.isna(r.bg_frac) else pc(r.bg_frac)} | {'' if pd.isna(r.p_depletion) else f'{r.p_depletion:.3f}'} | "
                 f"{'' if pd.isna(pv) else f'{pv:.3g}'} |")
    g = sev.set_index(["model", "test", "qc"])
    L += ["", f"**Agreement.** {int(g.loc[('raw','G2_lower_count','all'),'value'])} CS-lower genes are lower in Group 2 "
          f"({g.loc[('raw','G2_lower_count','all'),'genes']}). "
          f"{int(g.loc[('raw->cm_adjusted','G2_lower_persisting','all'),'value'])} of these remain lower after CM-content adjustment.",
          f"**Disagreement.** {int(g.loc[('raw','G2_higher_count','all'),'value'])} CS-lower genes are higher in Group 2 "
          f"({pc(g.loc[('raw','G2_higher_count','all'),'frac'])} vs background {pc(g.loc[('raw','G2_higher_count','all'),'bg_frac'])}, "
          f"enrichment p = {g.loc[('raw','G2_higher_count','all'),'p_enrichment']:.3f}): "
          f"{g.loc[('raw','G2_higher_count','all'),'genes']}. "
          f"{int(g.loc[('raw->cm_adjusted','G2_higher_persisting','all'),'value'])} remain higher after adjustment. "
          f"Across all {int(g.loc[('raw','spearman_all_overlap','all'),'n'])} overlapping genes, spatial log2 FC and the G2−G1 effect "
          f"correlate at ρ = {g.loc[('raw','spearman_all_overlap','all'),'value']:+.3f}.",
          "The proteome is bulk tissue that includes fibrotic regions, whereas the spatial spots are CM-dominant by construction; "
          "per the handoff, the ECM half is not interpreted as biology.", ""]
    return L


L = ["# ESC HF short report: results summary", "",
     "Generated by `scripts/eschf/a8_summary.py` from the ESC HF tables. Every value here is also checked by",
     "`scripts/manuscript_results.py` (section 4 of `docs/DISCREPANCIES.md`). Interpretation is limited to what each",
     "item requires.", ""]

# 0. revision: new primary ---------------------------------------------------
cmp_ = rdp(TOP, "a0_vs_old_primary.tsv").set_index("gene_set")
U = rdp(TOP, "a0_universe.tsv")
pm = rdp(V2, "section_patient_map.tsv")
ret = pm[pm.retained == "yes"]
excl = pm[(pm.arm == "CS") & (pm.retained != "yes")].section.tolist()
c2 = rdp(V2, "a1_patient_counts.tsv").set_index("arm")
L += ["## 0. Primary contrast (revision): all cardiomyocyte-dominant spots", "",
      f"The primary contrast is CS vs comparators over all CM-dominant spots, with no exclusion by distance to immune-dominant "
      f"spots: {int((ret.arm == 'CS').sum())} CS sections from {ret[ret.arm == 'CS'].patient.nunique()} patients "
      f"({c2.loc['CS','spots_retained']} spots) and {int((ret.arm == 'COMPARATOR').sum())} comparator sections from "
      f"{ret[ret.arm == 'COMPARATOR'].patient.nunique()} patients ({c2.loc['COMPARATOR','spots_retained']} spots); "
      f"{', '.join(excl)} has no CM-dominant spot and is not in the fit. "
      "The earlier fit on CM-dominant spots ≥0.5 mm from an immune-dominant spot is reported as a sensitivity analysis (Part II).", "",
      f"Gene universe: the committed rule (detection fraction ≥0.25 in CS or in comparators), applied to the new spot set: "
      f"{int(U.in_new_universe.sum())} genes (earlier fit {int(U.in_old_universe.sum())}; "
      f"{int((U.in_new_universe & U.in_old_universe).sum())} in both).", "",
      "| Gene set | New hits (down/up) | Earlier hits (down/up) | In both | New only | Earlier only | Earlier hits outside new universe | Earlier hits same sign in new fit | Pearson r, logFC (shared genes) |",
      "|---|---|---|---|---|---|---|---|---|"]
for r in cmp_.itertuples():
    gs = r.Index
    L.append(f"| {gs} | {r.n_hits_new} ({r.n_hits_new_down}/{r.n_hits_new_up}) | {r.n_hits_old} ({r.n_hits_old_down}/{r.n_hits_old_up}) | "
             f"{r.n_overlap} | {r.n_new_only} | {r.n_old_only} | {r.n_old_hits_not_in_new_universe} | "
             f"{r.old_hits_same_sign_in_new}/{r.n_old_hits_in_new_universe} | {r.pearson_r_logFC_shared:.3f} |")
prog = [l.strip() for l in open(f"{V2}/geneset_primary_program_qcpass_down.txt") if l.strip() and not l.startswith("#")]
nconc2 = len(rdp(V2, "a4_snrna_nf_concordant18.tsv")); nconc1 = len(rdp(TOP, "a4_snrna_nf_concordant18.tsv"))
c_new = set(rdp(V2, "a4_snrna_nf_concordant18.tsv").gene); c_old = set(rdp(TOP, "a4_snrna_nf_concordant18.tsv").gene)
L += ["", f"Program score genes (QC-passing, lower-in-CS hits of the new fit): {len(prog)}. Concordant genes "
      f"(v2 Panel E definition, section-stable core): {nconc2} ({len(c_new & c_old)} shared with the {nconc1} of the earlier fit).", "",
      "**Circularity.** The program genes are hits of the primary fit, which now includes the spots nearest immune-dominant "
      "spots, so the level difference of the program score in the near bins is not an out-of-sample test. The OXPHOS "
      "level check (below) uses a gene set defined independently of the fit.", ""]

sl = rdp(TOP, "a2_perpatient_slopes.tsv").set_index("outcome")
L += ["### 0.1 Per-patient slopes on distance (pre-specified, handoff s.3)", "",
      "For each CS patient with CM-dominant spots in at least two populated bins: OLS on that patient's spots, "
      "outcome ~ distance + CM-structural score + z(log10 UMI) (+ section indicator for patients with >1 section); "
      "the positive control (CM-structural score) is fitted without the structural covariate. Mean of per-patient slopes, "
      "95% t-interval.", "",
      "| Outcome | Patients | Mean slope per mm | 95% CI | Negative / positive | Gap share across max distance |", "|---|---|---|---|---|---|"]
for r in sl.itertuples():
    k = r.Index
    L.append(f"| {k} | {r.n_patients} | {r.mean_slope_per_mm:+.4f} | {r.ci_low:+.4f} to {r.ci_high:+.4f} | {r.n_negative}/{r.n_positive} | "
             f"{'' if pd.isna(r.gap_share_across_max_distance) else f3(r.gap_share_across_max_distance)} |")
L += ["", f"Gap share = |mean slope × {f2(sl.max_distance_mm.iloc[0])} mm| / (comparator mean of patient means − CS mean of "
      f"patient means). Maximum distance from any CS CM-dominant spot to the nearest immune-dominant spot: "
      f"{f2(sl.max_distance_mm.iloc[0])} mm; the 2–3 mm, >3 mm and no-immune-dominant-spot bins are empty.", ""]

ox = rdp(TOP, "a2_oxphos_level.tsv")
L += ["### 0.2 OXPHOS level check by bin (KEGG hsa00190; handoff s.4)", "",
      "| Bin | Spots | Patients | CS mean of patient means (95% CI) | Comparator mean (range) | CS mean below comparator range |",
      "|---|---|---|---|---|---|"]
for r in ox.itertuples():
    L.append(f"| {r.bin} | {r.n_spots} | {r.n_patients} | {f3(r.cs_mean_of_patient_means)} ({f3(r.cs_ci_low)}–{f3(r.cs_ci_high)}) | "
             f"{f3(r.comparator_mean)} ({f3(r.comparator_min)}–{f3(r.comparator_max)}) | {r.cs_mean_below_comparator_min} |")
L += ["", f"In {int((~ox.cs_mean_below_comparator_min.astype(bool)).sum())} of {len(ox)} populated bins the CS mean lies within "
      "or above the comparator range. Whether this table enters the figure, legend or text is Rami's decision (handoff s.4).", ""]

dd = rdp(V2, "a3a4_direction_descriptive.tsv")
L += ["### 0.3 Direction agreement, lower-in-CS hits (descriptive; handoff s.6)", "",
      "Each comparison is reported as the fraction in the stated direction beside the background fraction; p-values are given "
      "for all four comparisons.", "",
      "| Comparison | Gene set | Stated direction | Testable | In direction | Background | p (method) |", "|---|---|---|---|---|---|---|"]
for r in dd.itertuples():
    L.append(f"| {r.comparison} | {r.gene_set} | {r.stated_direction} | {r.n_testable} | {r.n_in_direction} ({pc(r.frac_in_direction)}) | "
             f"{pc(r.background_frac)} | {r.p_value:.3g} ({r.p_method}) |")
pv = rdp(TOP, "a7_panel_values.tsv")
def _pvD(q):
    v = pv[(pv.panel == "D") & (pv.quantity == q)].value
    assert len(v) == 1, q
    return float(v.iloc[0])
assert _pvD("max_abs_snrna_vs_nf_logFC") <= 3, "NF column exceeds the ±3 screen"
oor = pv[(pv.panel == "C") & pv.quantity.str.startswith("out_of_range|") & pv.quantity.str.endswith("|snrna_lfc")]
assert len(oor) == int(pv[(pv.panel == "C") & (pv.quantity == "n_out_of_range")].value.iloc[0]), "triangle list != off-scale count"
L += ["", f"Figure 1C: y-axis limited to ±3; {len(oor)} genes beyond it are drawn as outward triangles: "
      + ", ".join(f"*{q.split('|')[1]}* ({v:+.1f})" for q, v in zip(oor.quantity, oor.value)) + ".", ""]

# Part I: primary --------------------------------------------------------------
L += ["## Part I. Primary fit: all CM-dominant spots", ""]
L += fit_sections(V2, f"{TOP}/a0_primary_refit_summary.tsv", rdp(V2, "q1_hits_annotated.tsv"), "primary", 1,
                  f"{V2}/geneset_primary_program_qcpass_down.txt")
a6 = rdp(V2, "a6_tissue_fibrosis.tsv")
L += ["### 1.6 A6: tissue-level fibrosis, CS vs comparators (all spots, patient means; independent of the fit)", "",
      "| Variable | CS mean (median), n | Comparator mean (median), n | CS − comparator (95% CI) | Welch p | Mann–Whitney p |",
      "|---|---|---|---|---|---|"]
for _, r in a6.iterrows():
    L.append(f"| {r.variable} | {f3(r.mean_CS)} ({f3(r.median_CS)}), {r.n_CS_patients} | {f3(r.mean_comparator)} "
             f"({f3(r.median_comparator)}), {r.n_comparator_patients} | {r.diff_CS_minus_comparator:+.3f} "
             f"({r.ci_low:+.3f} to {r.ci_high:+.3f}) | {r.welch_p:.3f} | {r.mannwhitney_p:.3f} |")
L += ["", "No Kuppe non-failing sections are in the export, so no normal reference was computed.", ""]

# Part II: sensitivity ---------------------------------------------------------
L += ["## Part II. Sensitivity analysis: CM-dominant spots ≥0.5 mm from an immune-dominant spot (earlier primary)", ""]
L += fit_sections(TOP, f"{TOP}/a1_patient_refit_summary.tsv", rdp(TOP, "q1_hits_annotated.tsv"), "sensitivity, ≥0.5 mm", 2,
                  "metadata/geneset_primary_program_qcpass_down.txt")

# Part III: specifications not run as written ---------------------------------
h2 = rdp(V2, "q1_hits_annotated.tsv"); h2 = h2[h2["adj.P.Val"] < FDR]
m2 = rdp(V2, "a2_distance_model.tsv")
L += ["## 3. Specifications that could not be run as written", "",
      f"1. A2 bins 2–3 mm, >3 mm and no-immune-dominant-spot: empty (maximum distance {f2(sl.max_distance_mm.iloc[0])} mm; every "
      "CS section contains immune-dominant spots). Panel B shows populated bins only.",
      "2. A2/A6 normal-myocardium reference: the per-spot export contains CS and comparator sections only (no Kuppe normals).",
      f"3. A2 spot-level mixed models: {int(truth(m2.singular).sum())} of {len(m2)} fits are singular in the primary. They are "
      "reported as fitted, labelled as spot-level with anti-conservative intervals; the per-patient slopes are the pre-specified summary.",
      "4. A5 sarcomeric content score: *ACTC1* is not in the Hamed proteome. The score uses the remaining seven markers, with no substitute.",
      "5. A5b disease duration: not in the Neyazi supplementary tables S1/S2.",
      "6. Handoff label \"absent annotation\": the affected genes are those outside the locus-QC universe "
      "(`in_qc_universe == False`), not genes with `flag_missing_annotation` "
      f"(set for {int(truth(h2.flag_missing_annotation).sum()) if 'flag_missing_annotation' in h2 else 'none'} of the {len(h2)} primary hits).",
      "7. The IDH2 genotype-strata claim is dropped, as instructed.",
      "8. Universe rule: handoff s.0 describes detection in CS; the committed rule (and the one applied) is detection ≥0.25 in "
      "CS or in comparators (`01_limma_q1q2.R`).",
      f"9. {', '.join(excl)} has no CM-dominant spot, so the primary fit has {int((ret.arm == 'CS').sum())} CS sections rather than "
      f"{int((pm.arm == 'CS').sum())}.",
      "10. A12: *KCTD14* (NDUFC2 locus partner) is not in the snRNA count matrices; its row in the suspect-loci table is empty.",
      "11. A11: apex sampling cannot be recorded; Table S2 has no apex category. One S1 CS subject (SNU-1) has no snRNA sample "
      "in the GEO download, and its two S2 breakdowns disagree (0 by region, 2 by histology).",
      f"12. Figure 1D (v5) uses the {int(_pvD('n_concordant_genes'))} concordant genes of the section-stable core after screen 3a "
      f"(user decision, replacing the earlier all-hits set); {int(_pvD('n_core_removed_screen_3a'))} core gene(s) failing 3a removed. Panel C draws all screen-3a failures as "
      "open grey triangles. "
      f"The ±3 log2 screen was also applied to the snRNA vs NF column of Figure 1D: no shown gene exceeds it (largest "
      f"|log2FC| {_pvD('max_abs_snrna_vs_nf_logFC'):.2f}; lowest mean expression {_pvD('min_snrna_vs_nf_AveExpr'):.2f} log2 CPM), "
      "so it removes no further gene. "
      "Screen 3b is reported and not applied.", ""]

# Part IV: missing genes (primary) ---------------------------------------------
L += ["## 4. Genes missing from a dataset, by analysis (primary fit)", ""]
pool = rdp(V2, "snrna_pooled_full.tsv").set_index("gene")
miss_sn = sorted(set(h2.gene) - set(pool.index[np.isfinite(pool.lfc_pool)]))
a2m = pd.read_csv("data/eschf/primary_v2/a2_missing_genes.tsv", sep="\t")
a5m = rdp(V2, "a5_missing_genes.tsv"); a5bm = rdp(V2, "a5b_missing_genes.tsv")
L += [f"- A3 / A5b cross-platform (spatial hits not testable in pooled snRNA, {len(miss_sn)}): {', '.join(miss_sn) or 'none'}.",
      "- A2 scores (absent from the Visium gene space): "
      + ("; ".join(f"{r.score}: {r.gene}" for _, r in a2m.iterrows()) if len(a2m) else "none") + ".",
      f"- A5b snRNA program score (absent from the snRNA count matrix, {len(a5bm)}): {', '.join(a5bm.gene) or 'none'}."]
for (an, gs), d in a5m.groupby(["analysis", "gene_set"], sort=False):
    L.append(f"- {an} ({gs}, {len(d)}): {', '.join(sorted(d.gene))}.")
ex = rdp(V2, "a3_excluded_genes.tsv")
L += [f"- A3 QC-pass subset: {len(ex)} hits excluded, by reason: "
      + ", ".join(f"{k} {v}" for k, v in ex.exclusion_reason.value_counts().items()) + ".", ""]

# Part V (written last): revision v3 (HANDOFF_figure_v3_ESCHF.md) ------------------------------
L += ["## 5. Revision v3: snRNA reliability screen, pathway levels, cohort accounting, GPNMB", ""]
S9 = rdp(TOP, "a9_snrna_screen.tsv").set_index("gene")
cnt = rdp(TOP, "a9_screen_counts.tsv").set_index("set")
L += ["### 5.1 A9 snRNA reliability screen (pre-specified)", "",
      "3a (applied to Figure 1C/D): a spatial hit is unreliable in snRNA if |pooled snRNA log2FC| > 3. "
      "3b (reported, not applied): source-specific near-absence, i.e. some snRNA study has per-study mean log2 CPM < 0 while "
      "the mean of the other four studies' means is ≥ 2 higher. The threshold is a pre-stated screen.", "",
      "| Set | Genes | Fail 3a | 3a genes | Flagged 3b | 3a and 3b | 3a only | 3b only |", "|---|---|---|---|---|---|---|---|"]
for k in ["all_hits_testable", "qc_pass_hits_testable", "concordant_allhits", "concordant_core_v2"]:
    r = cnt.loc[k]
    L.append(f"| {k} | {int(r.n)} | {int(r.n_fail_3a)} | {', '.join(f'*{g}*' for g in str(r.genes_fail_3a).split(','))} | "
             f"{int(r.n_flag_3b)} | {int(r.n_3a_and_3b)} | {int(r.n_3a_only)} | {int(r.n_3b_only)} |")
L += ["", "Published Chaffin 2022 cardiomyocyte contrasts vs non-failing, genes with |log2FC| > 3: "
      + "; ".join(f"{t} {int(cnt.loc[f'chaffin_{t}_vs_NF_published', 'n_fail_3a'])} of {int(cnt.loc[f'chaffin_{t}_vs_NF_published', 'n'])} "
                  f"({cnt.loc[f'chaffin_{t}_vs_NF_published', 'genes_fail_3a']})" for t in ("DCM", "HCM")) + ".", ""]
sa = rdp(TOP, "a9_source_absence.tsv")
studs = [c for c in sa.columns if c.startswith("mean_log2cpm_")]
L += ["3b-flagged hits (and the FPGT reference gene), per-study mean log2 CPM:", "",
      "| Gene | Direction | Spatial log2FC | Pooled snRNA log2FC | Fails 3a | Absent in | " + " | ".join(c.replace("mean_log2cpm_", "") for c in studs) + " |",
      "|---|---|---|---|---|---|" + "---|" * len(studs)]
for r in sa.itertuples():
    sl_ = "" if pd.isna(r.spatial_logFC) else f"{r.spatial_logFC:+.2f}"
    lp_ = "" if pd.isna(r.lfc_pool) else f"{r.lfc_pool:+.2f}"
    L.append(f"| *{r.gene}* | {'' if pd.isna(r.direction) else r.direction} | {sl_} | {lp_} | "
             f"{'' if pd.isna(r.fail_3a) else r.fail_3a} | {r.absent_study} | " + " | ".join(f2(getattr(r, c)) for c in studs) + " |")
lo = rdp(TOP, "a9_suspect_loci.tsv")
L += ["", "Suspect loci and partners (per-study mean log2 CPM; KCTD14 is not in the snRNA count matrices):", "",
      "| Locus | Role | " + " | ".join(c.replace("mean_log2cpm_", "") for c in studs) + " | Pooled snRNA log2FC |",
      "|---|---|" + "---|" * len(studs) + "---|"]
for _, r in lo.iterrows():
    L.append(f"| {r.gene} | {r.role} | " + " | ".join("" if pd.isna(r[c]) else f2(r[c]) for c in studs)
             + f" | {'' if pd.isna(r.lfc_pool) else f'{r.lfc_pool:+.2f}'} |")
cc = rdp(TOP, "a9_concordance_screened.tsv")
L += ["", "Sign agreement with the pooled snRNA contrast (3c; background from the expression-decile-matched label-permutation "
      "null, 2000 permutations; screened = 3a failures removed from hits and null pool):", "",
      "| Screen | Direction | Gene set | Testable | Same sign | Background | Same sign, snRNA FDR <5% |", "|---|---|---|---|---|---|---|"]
for r in cc.itertuples():
    L.append(f"| {r.screen} | {r.direction} | {r.gene_set} | {r.n_testable} | {r.n_same_sign} ({pc(r.frac_same_sign)}) | "
             f"{pc(r.background_frac)} | {r.n_concordant_fdr5} |")
ca = S9.index[truth(S9.in_concordant_allhits_screened)]; cv = S9.index[truth(S9.in_concordant_core_v2)]
ae = S9.loc["ENG"]
L += ["", f"Concordant genes, v3 Figure 1D definition (all spatial hits, same sign, snRNA FDR <5%, passing 3a): {len(ca)}. "
      f"The v2 panel used the section-stable core ({len(cv)} genes); {len(set(ca) & set(cv))} are in both. "
      f"*ENG*: same sign {bool(ae.same_sign)}, snRNA FDR {ae.fdr_pool:.2g}, pooled log2FC {ae.lfc_pool:+.2f}; it is not in the "
      f"section-stable core (in_concordant_core_v2 = {bool(truth(pd.Series([ae.in_concordant_core_v2])).iloc[0])}), which is why it "
      "was absent from the v2 concordant panel, and it fails 3a.", ""]
k12 = rdp(TOP, "a12_key_genes.tsv")
L += ["### 5.2 A12 key genes", "", "| Gene | Spatial log2FC | Pooled snRNA log2FC (FDR) | DCM | ARVC | Chaffin DCM | Chaffin HCM | QC pass | 3a fail | 3b flag | Concordant (all hits, screened) |",
      "|---|---|---|---|---|---|---|---|---|---|---|"]
fm = lambda x: "" if pd.isna(x) else f"{x:+.2f}"
for r in k12.itertuples():
    L.append(f"| *{r.gene}* | {fm(r.spatial_logFC)} | {fm(r.lfc_pool)} ({r.fdr_pool:.2g}) | {fm(r.lfc_DCM)} | {fm(r.lfc_ARVC)} | "
             f"{fm(r.chaffin_DCM_logFC)} | {fm(r.chaffin_HCM_logFC)} | {r.qc_pass} | {r.fail_3a} | {r.flag_3b} | {r.in_concordant_allhits_screened} |")
a10 = rdp(TOP, "a10_pathway_levels.tsv"); m10 = rdp(TOP, "a10_pathway_members.tsv")
L += ["", "### 5.3 A10 a priori pathway levels (KEGG; gene sets defined independently of the fit)", "",
      "Spot level: mean log1p CP10K over in-universe members on depth-matched counts, CS mean of patient means per bin vs the "
      "comparator range of per-patient means. Patient level: mean voom-form log2 CPM over members on patient pseudobulk, CS − "
      "comparator, Welch t-test (95% CI).", "",
      "| Pathway | Level / bin | Genes | CS | Comparator mean (range) | CS below comparator min | CS − comparator (95% CI) | Welch p |",
      "|---|---|---|---|---|---|---|---|"]
for r in a10.itertuples():
    pat = r.level == "patient_pseudobulk"
    cs_ = f3(r.cs_mean_of_patient_means) + ("" if pd.isna(r.cs_ci_low) else f" ({f3(r.cs_ci_low)}–{f3(r.cs_ci_high)})")
    L.append(f"| {r.set} | {'patient pseudobulk' if pat else r.bin} | {r.n_genes} | {cs_} | {f3(r.comparator_mean)} "
             f"({f3(r.comparator_min)}–{f3(r.comparator_max)}) | {'' if pat else r.cs_mean_below_comparator_min} | "
             + (f"{r.diff_cs_minus_comparator:+.3f} ({r.diff_ci_low:+.3f} to {r.diff_ci_high:+.3f}) | {r.welch_p:.3g} |" if pat else " | |"))
ov = m10[truth(m10.in_universe)].groupby("set").agg(n=("gene", "size"), hits=("primary_hit", "sum"),
                                                     down=("hit_direction", lambda v: int((v == "down").sum())))
L += ["", "Overlap with primary-fit hits (in-universe members): " + "; ".join(
      f"{k} {int(r.hits)} of {int(r.n)} ({int(r.down)} lower in CS)" for k, r in ov.iterrows())
      + ". These pathway scores therefore share genes with the hit list.", ""]
a11 = rdp(TOP, "a11_snrna_cohort_accounting.tsv")
an = truth(a11.analysed)
L += ["### 5.4 A11 snRNA cohort accounting (Neyazi Tables S1/S2 vs the analysed cohort)", "",
      f"S1 CS subjects: {len(a11)}; analysed: {int(an.sum())}; specimen count equals the S2 region count for "
      f"{int(truth(a11[an].specimens_match_s2_region).sum())} of {int(an.sum())} analysed patients. Analysed by tissue source: "
      + ", ".join(f"{k} {v}" for k, v in a11[an]["Tissue Source"].value_counts().items()) + ". "
      "Apex is not a Table S2 category, so apical sampling is not recorded.", "",
      "| Subject | Tissue source | Analysed | S2 regions sampled | Specimens analysed | S2 by region | S2 by histology | Reason not analysed |",
      "|---|---|---|---|---|---|---|---|"]
for r in a11.itertuples():
    L.append(f"| {r.Subject} | {getattr(r, '_3')} | {r.analysed} | {r.s2_regions_sampled} | "
             f"{'' if pd.isna(r.n_specimens) else int(r.n_specimens)} | {r.s2_n_by_region} | {r.s2_n_by_histology} | "
             f"{'' if pd.isna(r.missing_reason) else r.missing_reason} |")
# A13 (GPNMB source check) removed: it depended on the Kuppe normal sections, which are no longer in the export.

# Part VI (written last): A14 pathway scores with fit hits removed ---------------
R14 = rdp(TOP, "a14_pathway_hits_removed.tsv"); N14 = rdp(TOP, "a14_pathway_null.tsv").set_index("set")
B14 = rdp(TOP, "a14_pathway_robustness.tsv")
PNAME = {"tca_cycle_hsa00020": "TCA cycle (hsa00020)", "fatty_acid_degradation_hsa00071": "Fatty acid degradation (hsa00071)",
         "oxphos_hsa00190": "Oxidative phosphorylation (hsa00190)"}
L += ["## 6. A14: a priori pathway scores with fit hits removed (HANDOFF_pathway_hits_removed.md)", "",
      "Rule (pre-stated): the hits-removed (primary) difference holds if its Welch 95% interval excludes zero and its "
      "percentile in the expression-matched random-set null is below 2.5 or above 97.5. The rule does not specify a "
      "direction; direction is stated separately.", ""]
for k, nm in PNAME.items():
    r = R14[(R14.set == k) & (R14.version == "minus_primary")].iloc[0]; f = R14[(R14.set == k) & (R14.version == "full")].iloc[0]
    n = N14.loc[k]; ex0 = r.diff_ci_low > 0 or r.diff_ci_high < 0; tail = n.percentile < 2.5 or n.percentile > 97.5
    L.append(f"- **{nm}**: with primary-fit hits removed ({r.n_genes} genes), CS − comparator {r.diff_cs_minus_comparator:+.3f} "
             f"log2 CPM (95% CI {r.diff_ci_low:+.3f} to {r.diff_ci_high:+.3f}), null percentile {n.percentile:.1f}; the interval "
             f"{'excludes' if ex0 else 'includes'} zero and the percentile is {'beyond' if tail else 'within'} 2.5–97.5, so the "
             f"difference {'holds' if ex0 and tail else 'does not hold'} by the stated rule. CS is {'higher' if r.diff_cs_minus_comparator > 0 else 'lower'} "
             f"than comparators in this version; the full set ({f.n_genes} genes, A10) gave {f.diff_cs_minus_comparator:+.3f} "
             f"({f.diff_ci_low:+.3f} to {f.diff_ci_high:+.3f}).")
L += ["", "### 6.1 Gene-set versions and patient-level difference (8 CS vs 6 comparator patients; mean voom-form log2 CPM)", "",
      "| Pathway | Version | Genes | Removed | CS mean | Comparator mean (range) | CS − comparator (95% CI) | Welch p | Too small |",
      "|---|---|---|---|---|---|---|---|---|"]
for r in R14.itertuples():
    rem = "" if pd.isna(r.removed_genes) else ", ".join(f"*{g}*" for g in r.removed_genes.split(","))
    L.append(f"| {PNAME[r.set]} | {r.version} | {r.n_genes} | {rem} | {f3(r.cs_mean)} | {f3(r.comparator_mean)} "
             f"({f3(r.comparator_min)}–{f3(r.comparator_max)}) | {r.diff_cs_minus_comparator:+.3f} ({r.diff_ci_low:+.3f} to "
             f"{r.diff_ci_high:+.3f}) | {r.welch_p:.3g} | {'yes' if truth(pd.Series([r.too_small_to_interpret])).iloc[0] else 'no'} |")
L += ["", "### 6.2 Expression-matched random-set null (hits-removed primary version)", "",
      "Pool: tested-universe genes minus all primary-fit hits minus all members of the pathway tested. Expression = mean "
      "voom-form log2 CPM over the 14 patients; deciles over the 982-gene universe; each draw matches the members' decile "
      "counts. 10,000 draws, seed 0.", "",
      "| Pathway | Genes | Pool | Observed | Null mean (SD) | Null 2.5–97.5% | Percentile | Two-sided p |", "|---|---|---|---|---|---|---|---|"]
for k, n in N14.iterrows():
    L.append(f"| {PNAME[k]} | {n.n_genes} | {n.n_pool} | {n.observed_diff:+.3f} | {n.null_mean:+.3f} ({f3(n.null_sd)}) | "
             f"{n.null_q025:+.3f} to {n.null_q975:+.3f} | {n.percentile:.1f} | {n.p_two_sided:.3g} |")
L += ["", "### 6.3 Robustness (hits-removed primary version)", "",
      "| Pathway | Leave-one-gene-out range | Most influential gene (change) | Members lower / higher in CS | Per-patient slope per mm (95% CI; negative/positive) | Slope t p |",
      "|---|---|---|---|---|---|"]
for k, nm in PNAME.items():
    b = B14[B14.set == k]
    lo = b[b.analysis == "leave_one_gene_out_summary"].iloc[0]; dc = b[b.analysis == "per_gene_direction_counts"].iloc[0]
    ss = b[b.analysis == "slope_summary"].iloc[0]
    L.append(f"| {nm} | {lo.value:+.3f} to {lo.value2:+.3f} | *{lo.gene}* ({lo.change:+.3f}) | {int(dc.value)} / {int(dc.value2)} | "
             f"{ss.value:+.4f} ({ss.ci_low:+.4f} to {ss.ci_high:+.4f}; {int(ss.n_negative)}/{int(ss.n_positive)}) | {ss.t_p:.3g} |")
L += ["", "Spot level by distance (mean log1p CP10K over members, T = 2000; CS mean of patient means with 95% CI vs comparator patient means):", "",
      "| Pathway | Bin | CS patients | CS mean (95% CI) | Comparator mean (range) | CS mean below comparator min |", "|---|---|---|---|---|---|"]
for r in B14[B14.analysis == "distance_bin"].itertuples():
    L.append(f"| {PNAME[r.set]} | {r.item} | {int(r.n_patients)} | {f3(r.value)} ({f3(r.ci_low)}–{f3(r.ci_high)}) | "
             f"{f3(r.comparator_mean)} ({f3(r.comparator_min)}–{f3(r.comparator_max)}) | {r.cs_mean_below_comparator_min} |")
L += ["", "Per-member CS − comparator values and leave-one-gene-out values per gene: `results/eschf/a14_pathway_robustness.tsv`; "
      "every null draw's genes: `a14_pathway_null_draws.tsv.gz`.", ""]

open(f"{TOP}/SUMMARY.md", "w").write("\n".join(L))
print(len(L), "lines")
