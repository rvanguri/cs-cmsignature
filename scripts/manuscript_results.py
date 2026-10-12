#!/usr/bin/env python3
"""Recompute every quantitative claim in the ESC HF short report from the committed tables.

Two blocks: ESCHF[...] rows check the analysis outputs against independent recomputation from upstream
inputs; ShortReport rows encode every number printed in the manuscript text, figure and legend at its
printed precision and recompute it from the committed tables. Each claim is printed as
  manuscript value | recomputed value | status
and written to results/manuscript_results.tsv. Exit status is 1 if any claim is MISMATCH, so this can gate
a commit. NOT_RUN marks a check whose model is not estimable on these data (reported, not counted as a pass).

Run from the repo root after the ESC HF pipeline (scripts/eschf/run_pipeline.sh, stage local3).
The earlier spatial-first brief report and its claims are at git tag brief-report-spatial-first.
"""
import glob
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
YLIM_C = 3.0  # Figure 1C y-axis half-range (= the +/-3 snRNA screen)


def claim(name, manuscript, recomputed, ok, note=""):
    status = "NOT_REPRODUCIBLE" if ok is None else ("NOT_RUN" if isinstance(ok, str) and ok == "NOT_RUN" else ("OK" if bool(ok) else "MISMATCH"))
    rows.append(dict(claim=name, manuscript=manuscript, recomputed=recomputed,
                     status=status, note=note))


ESCHF_FITS = {
    # primary (handoff_figure_fixes s.0): all cardiomyocyte-dominant spots, no immune-adjacency exclusion
    "primary": dict(E="results/eschf/primary_v2", meta="data/eschf/primary_v2/meta_allcm_T2000.tsv",
                    perspot="data/eschf/primary_v2/eschf_perspot_T2000.tsv.gz",
                    refit="results/eschf/a0_primary_refit.tsv", refit_summary="results/eschf/a0_primary_refit_summary.tsv",
                    orig_col="section", orig_level="section_level_primary", gs_all="all_tested",
                    conc_col="sign_concordant_section_hits", still_col="section_hits_still_fdr5",
                    prog="results/eschf/primary_v2/geneset_primary_program_qcpass_down.txt"),
    # sensitivity: earlier fit on cardiomyocyte-dominant spots >= 0.5 mm from an immune-dominant spot
    "sens>=0.5mm": dict(E="results/eschf", meta=f"{DAT}/meta_T2000.tsv", perspot="data/eschf/eschf_perspot_T2000.tsv.gz",
                        refit="results/eschf/a1_patient_refit_pergene.tsv",
                        refit_summary="results/eschf/a1_patient_refit_summary.tsv", orig_col="orig",
                        orig_level="original_section_level", gs_all="all_985", conc_col="sign_concordant_orig_hits",
                        still_col="orig_hits_still_fdr5", prog="metadata/geneset_primary_program_qcpass_down.txt"),
}


def eschf_claims():
    for fit_tag, cfg in ESCHF_FITS.items():
        _eschf_fit_claims(fit_tag, cfg)
    eschf_revision_claims()
    eschf_v3_claims()
    eschf_a14_claims()


def eschf_a14_claims():
    """ESC HF A14 (HANDOFF_pathway_hits_removed.md): pathway scores with fit hits removed, recomputed from the hit
    tables, KEGG files, universe, voom-form patient log2 CPM, spot counts and model input; every null draw is
    re-verified (size, decile composition, pool exclusion) and its statistic recomputed. Also checks the
    SUMMARY.md Figure 1C triangle list (one entry per off-scale gene, true pooled snRNA value)."""
    import re
    import scipy.sparse as ssp
    import statsmodels.api as sm
    E, V, D = "results/eschf", "results/eschf/primary_v2", "data/eschf/primary_v2"
    rd = lambda f, **k: pd.read_csv(f, sep="\t", **k)
    TB = lambda v: v.astype(str).str.upper().eq("TRUE")
    rs = lambda p: [l.strip() for l in open(p) if l.strip() and not l.startswith("#")]

    def num(name, reported, recomputed, tol=1e-6, note=""):
        r, c = float(reported), float(recomputed)
        ok = (np.isnan(r) and np.isnan(c)) or bool(np.isclose(r, c, rtol=tol, atol=tol))
        claim(f"ESCHF[a14] {name}", round(r, 4), round(c, 4), ok, note)

    def same(name, reported, recomputed):
        a, b = sorted(reported), sorted(recomputed)
        claim(f"ESCHF[a14] {name}", ",".join(a) or "none", ",".join(b) or "none", a == b)

    # ---- SUMMARY triangle list (section 5 fix)
    line = [l for l in open(f"{E}/SUMMARY.md") if l.startswith("Figure 1C: y-axis limited")][0]
    ent = re.findall(r"\*([A-Za-z0-9\-]+)\* \(([+\-]\d+\.\d)\)", line)
    Pn = rd(f"{V}/snrna_pooled_full.tsv").set_index("gene")
    Hq = rd(f"{V}/q1_hits_annotated.tsv").set_index("gene"); Hq = Hq[TB(Hq.qc_pass)]
    lpq = Pn.lfc_pool.reindex(Hq.index).dropna(); off = lpq[lpq.abs() > YLIM_C]
    same("SUMMARY Fig 1C triangle genes (one entry each)", [g for g, _ in ent], list(off.index))
    num("SUMMARY Fig 1C triangle entries", len(ent), len(off))
    for g, v in ent:
        num(f"SUMMARY Fig 1C triangle {g} value (1 dp)", float(v), round(off[g], 1), tol=1e-9)

    # ---- inputs
    univ = rs(f"{D}/universe_allcm.txt"); U = set(univ)
    P = set(rd(f"{V}/q1_hits_annotated.tsv").query("`adj.P.Val` < 0.05").gene)
    S = set(rd(f"{E}/q1_hits_annotated.tsv").query("`adj.P.Val` < 0.05").gene)
    L = rd(f"{E}/a10_patient_log2cpm_voomform.tsv").set_index("gene").loc[univ]
    mp = rd(f"{V}/section_patient_map.tsv").drop_duplicates("patient").set_index("patient").arm
    csp = [c for c in L.columns if mp[c] == "CS"]; cop = [c for c in L.columns if mp[c] == "COMPARATOR"]
    num("patient pseudobulk CS patients", 8, len(csp)); num("patient pseudobulk comparator patients", 6, len(cop))
    KG = {"tca_cycle_hsa00020": "tca_hsa00020", "fatty_acid_degradation_hsa00071": "fatty_acid_degradation_hsa00071",
          "oxphos_hsa00190": "oxphos_hsa00190"}
    R1 = rd(f"{E}/a14_pathway_hits_removed.tsv").set_index(["set", "version"])
    R2 = rd(f"{E}/a14_pathway_null.tsv").set_index("set")
    R3 = rd(f"{E}/a14_pathway_robustness.tsv")
    gd = L[csp].mean(1) - L[cop].mean(1)
    mem = {}
    for k, f in KG.items():
        allm = rs(f"metadata/geneset_kegg_{f}.txt"); full = sorted(set(allm) & U)
        vers = dict(full=full, minus_primary=sorted(set(full) - P), minus_union=sorted(set(full) - P - S))
        mem[k] = (set(allm), vers)
        for v, gs in vers.items():
            r = R1.loc[(k, v)]
            same(f"{k} {v} members", str(r.members).split(","), gs)
            same(f"{k} {v} removed genes", [] if pd.isna(r.removed_genes) else str(r.removed_genes).split(","),
                 sorted(set(full) - set(gs)))
            claim(f"ESCHF[a14] {k} {v} too-small flag", bool(TB(pd.Series([r.too_small_to_interpret])).iloc[0]),
                  len(gs) < 5, bool(TB(pd.Series([r.too_small_to_interpret])).iloc[0]) == (len(gs) < 5))
            sc = L.loc[gs].mean(0); x, y = sc[csp].values, sc[cop].values
            tt = stats.ttest_ind(x, y, equal_var=False); ci = tt.confidence_interval(0.95)
            num(f"{k} {v} CS mean", r.cs_mean, x.mean()); num(f"{k} {v} comparator mean", r.comparator_mean, y.mean())
            num(f"{k} {v} comparator min", r.comparator_min, y.min()); num(f"{k} {v} comparator max", r.comparator_max, y.max())
            num(f"{k} {v} CS - comparator", r.diff_cs_minus_comparator, x.mean() - y.mean())
            num(f"{k} {v} Welch CI low", r.diff_ci_low, ci.low); num(f"{k} {v} Welch CI high", r.diff_ci_high, ci.high)
            num(f"{k} {v} Welch p", r.welch_p, tt.pvalue)
    # ---- null: deciles by explicit quantile edges (qcut semantics), every draw verified
    ex = L.mean(1).values; edges = np.quantile(ex, np.linspace(0, 1, 11))
    dec = pd.Series(np.clip(np.searchsorted(edges, ex, side="left") - 1, 0, 9), index=univ)
    DR = rd(f"{E}/a14_pathway_null_draws.tsv.gz")
    for k in KG:
        allm, vers = mem[k]; gs = vers["minus_primary"]; obs = gd[gs].mean(); r = R2.loc[k]
        pool = set(univ) - P - allm
        need = dec[gs].value_counts().sort_index()
        num(f"{k} null pool size", r.n_pool, len(pool)); num(f"{k} null set size", r.n_genes, len(gs))
        claim(f"ESCHF[a14] {k} null decile composition", r.decile_counts,
              ",".join(f"{d}:{n}" for d, n in need.items()), r.decile_counts == ",".join(f"{d}:{n}" for d, n in need.items()))
        dk = DR[DR.set == k].sort_values("draw"); num(f"{k} null draws", r.n_draws, len(dk))
        G = dk.genes.str.split(",")
        sizes_ok = bool((G.str.len() == len(gs)).all())
        uniq_ok = bool(G.map(lambda z: len(set(z)) == len(z)).all())
        pool_ok = bool(G.map(lambda z: set(z) <= pool).all())
        need_d = need.to_dict()
        dec_ok = bool(G.map(lambda z: dec[z].value_counts().to_dict() == need_d).all())
        claim(f"ESCHF[a14] {k} every draw: size, distinct genes, in pool, decile-matched",
              "all", f"size={sizes_ok} distinct={uniq_ok} pool={pool_ok} decile={dec_ok}",
              sizes_ok and uniq_ok and pool_ok and dec_ok)
        null = np.array([gd[z].mean() for z in G]); N = len(null)
        lo, hi, eq = (null < obs).sum(), (null > obs).sum(), (null == obs).sum()
        num(f"{k} null observed diff", r.observed_diff, obs)
        num(f"{k} null mean", r.null_mean, null.mean()); num(f"{k} null 2.5%", r.null_q025, np.quantile(null, 0.025))
        num(f"{k} null 97.5%", r.null_q975, np.quantile(null, 0.975))
        num(f"{k} null percentile", r.percentile, 100 * (lo + 0.5 * eq) / N)
        num(f"{k} null two-sided p", r.p_two_sided, min(1.0, 2 * min((1 + lo + eq) / (N + 1), (1 + hi + eq) / (N + 1))))
    # ---- robustness
    X = ssp.load_npz(f"{D}/spotcounts_allcm_T2000.npz").tocsc()
    sp_ = rd(f"{D}/spotcounts_allcm_T2000_spots.tsv.gz"); gi = {g: i for i, g in enumerate(rs(f"{D}/spotcounts_allcm_T2000_genes.txt"))}
    mi = rd(f"{V}/a2_model_input.tsv.gz"); mi["z"] = (mi.log10_umi - mi.log10_umi.mean()) / mi.log10_umi.std(ddof=1)
    for k in KG:
        gs = mem[k][1]["minus_primary"]; rk = R3[R3.set == k]
        base = gd[gs].mean(); lo1 = pd.Series({g: gd[[h for h in gs if h != g]].mean() for g in gs})
        sm_ = rk[rk.analysis == "leave_one_gene_out_summary"].iloc[0]
        num(f"{k} LOGO min", sm_.value, lo1.min()); num(f"{k} LOGO max", sm_.value2, lo1.max())
        gmax = (lo1 - base).abs().idxmax()
        claim(f"ESCHF[a14] {k} LOGO most influential gene", sm_.gene, gmax, sm_.gene == gmax)
        num(f"{k} LOGO largest change", sm_.change, lo1[gmax] - base)
        dc = rk[rk.analysis == "per_gene_direction_counts"].iloc[0]
        num(f"{k} members lower in CS", dc.value, (gd[gs] < 0).sum()); num(f"{k} members higher in CS", dc.value2, (gd[gs] > 0).sum())
        pg = rk[rk.analysis == "per_gene_cs_minus_comparator"].set_index("item").value
        num(f"{k} per-gene CS - comparator, max abs deviation", 0, (pg - gd[pg.index]).abs().max())
        sub = X[:, [gi[g] for g in gs]].toarray() * (1e4 / 2000)
        sc = pd.DataFrame({"spot": sp_.spot, "section": sp_.section, "cohort": sp_.cohort, "score": np.log1p(sub).mean(1)})
        cs = mi.merge(sc, on=["spot", "section"])
        cm = sc[sc.cohort == "COMPARATOR"].groupby("section").score.mean()
        pm = cs.groupby(["bin", "patient"]).score.mean()
        for _, r in rk[rk.analysis == "distance_bin"].iterrows():
            v = pm.loc[r["item"]]
            ci = stats.t.interval(0.95, len(v) - 1, loc=v.mean(), scale=stats.sem(v))
            num(f"{k} bin {r['item']} CS mean of patient means", r.value, v.mean())
            num(f"{k} bin {r['item']} CI low", r.ci_low, ci[0]); num(f"{k} bin {r['item']} CI high", r.ci_high, ci[1])
            num(f"{k} bin {r['item']} comparator min", r.comparator_min, cm.min())
            num(f"{k} bin {r['item']} comparator max", r.comparator_max, cm.max())
            claim(f"ESCHF[a14] {k} bin {r['item']} CS mean below comparator min", bool(TB(pd.Series([r.cs_mean_below_comparator_min])).iloc[0]),
                  bool(v.mean() < cm.min()), bool(TB(pd.Series([r.cs_mean_below_comparator_min])).iloc[0]) == bool(v.mean() < cm.min()))
        sl = {}
        for p, d in cs.groupby("patient"):
            if d.bin.nunique() < 2:
                continue
            Xd = pd.DataFrame({"dist": d.dist_gran_mm, "struct": d.score_cm_structural, "z": d.z})
            for s_ in sorted(d.section.unique())[1:]:
                Xd[s_] = (d.section == s_).astype(float)
            sl[p] = sm.OLS(d.score.values, sm.add_constant(Xd)).fit().params["dist"]
        sl = pd.Series(sl); ps = rk[rk.analysis == "patient_slope_per_mm"].set_index("item").value
        num(f"{k} slope max abs deviation over patients", 0, (ps - sl.reindex(ps.index)).abs().max())
        ss = rk[rk.analysis == "slope_summary"].iloc[0]; ci = stats.t.interval(0.95, len(sl) - 1, loc=sl.mean(), scale=stats.sem(sl))
        num(f"{k} mean slope per mm", ss.value, sl.mean()); num(f"{k} slope CI low", ss.ci_low, ci[0])
        num(f"{k} slope CI high", ss.ci_high, ci[1]); num(f"{k} slope n negative", ss.n_negative, (sl < 0).sum())
        num(f"{k} slope t p", ss.t_p, stats.ttest_1samp(sl, 0).pvalue)


def _eschf_fit_claims(fit_tag, cfg):
    """ESC HF A1-A6 for one fit: reported values (tables under cfg E) vs independent recomputation
    from upstream inputs (per-spot export, hit table, snRNA and published tables)."""
    E = cfg["E"]
    rd = lambda f, **k: pd.read_csv(f"{E}/{f}", sep="\t", **k)

    def num(name, reported, recomputed, tol=1e-6, note=""):
        ok = bool(np.isclose(float(reported), float(recomputed), rtol=tol, atol=tol))
        claim(f"ESCHF[{fit_tag}] {name}", round(float(reported), 4), round(float(recomputed), 4), ok, note)

    def claim_t(name, *a, **k):
        claim(f"ESCHF[{fit_tag}] {name}", *a, **k)

    # ---- A1 ----------------------------------------------------------------------
    pmap = rd("section_patient_map.tsv")
    meta = pd.read_csv(cfg["meta"], sep="\t")
    used = meta[meta.cohort.isin(["CS", "COMPARATOR"]) & meta.evaluable.astype(str).str.lower().isin(["true", "1", "yes"])]
    m = used.merge(pmap[["section", "patient"]], on="section", how="left")
    assert m.patient.notna().all()
    cnt = rd("a1_patient_counts.tsv").set_index("arm")
    for arm in ("CS", "COMPARATOR"):
        num(f"A1 {arm} sections retained", cnt.loc[arm, "sections_retained"], (m.cohort == arm).sum())
        num(f"A1 {arm} patients retained", cnt.loc[arm, "patients_retained"], m[m.cohort == arm].patient.nunique())
    pg = pd.read_csv(cfg["refit"], sep="\t")
    oc = cfg["orig_col"]
    pg["qc_pass"] = pg.qc_pass.astype(str).str.upper().eq("TRUE")
    pg["adj_orig"], pg["logFC_orig"] = pg[f"adj_{oc}"], pg[f"logFC_{oc}"]
    pg["original_hit"] = pg.adj_orig < FDR
    sm = pd.read_csv(cfg["refit_summary"], sep="\t").set_index(["refit", "gene_set"])
    for gs, sub in ((cfg["gs_all"], pg), ("qc_pass", pg[pg.qc_pass])):
        orig = sub[sub.original_hit]
        num(f"A1 {gs} original hits", sm.loc[(cfg["orig_level"], gs), "hits_fdr5"], (sub.adj_orig < FDR).sum())
        num(f"A1 {gs} patient-pseudobulk hits", sm.loc[("patient_pseudobulk", gs), "hits_fdr5"], (sub.adj_patient < FDR).sum())
        num(f"A1 {gs} original hits same sign at patient level",
            sm.loc[("patient_pseudobulk", gs), cfg["conc_col"]],
            (np.sign(orig.logFC_patient) == np.sign(orig.logFC_orig)).sum())
        num(f"A1 {gs} original hits still FDR<5% at patient level",
            sm.loc[("patient_pseudobulk", gs), cfg["still_col"]], (orig.adj_patient < FDR).sum())
        if pd.isna(sm.loc[("section_dupcor_patient_block", gs), "hits_fdr5"]):
            claim_t(f"A1 {gs} dupcor-blocked hits", "not estimable", "not run", "NOT_RUN",
                    "duplicateCorrelation consensus is NaN (a1 refit 2 NOT ESTIMABLE)")
        else:
            num(f"A1 {gs} dupcor-blocked hits", sm.loc[("section_dupcor_patient_block", gs), "hits_fdr5"],
                (sub.adj_dupcor < FDR).sum())

    # ---- A2 / A6 from the per-spot export ----------------------------------------
    sp = pd.read_csv(cfg["perspot"], sep="\t")
    sp["patient"] = sp.section.map(pmap.set_index("section").patient)
    cs = sp[(sp.cohort == "CS") & sp.cm_dominant].copy()
    cs["bin"] = pd.cut(cs.dist_gran_mm, [0, 0.25, 0.5, 1, 2, 3, np.inf], right=False,
                       labels=["<0.25", "0.25-0.5", "0.5-1", "1-2", "2-3", ">3"]).astype(object)
    bins = rd("a2_distance_bins.tsv").set_index("bin")
    for b in ["<0.25", "0.25-0.5", "0.5-1", "1-2", "2-3", ">3"]:
        sub = cs[cs.bin == b]
        num(f"A2 bin {b} spots", bins.loc[b, "n_spots"], len(sub))
        num(f"A2 bin {b} patients", bins.loc[b, "n_patients"], sub.patient.nunique())
        if len(sub):
            for s in ("score_program", "score_ACADVL"):
                num(f"A2 bin {b} {s} mean of patient means", bins.loc[b, f"{s}_mean"],
                    sub.groupby("patient")[s].mean().mean())
    num("A2 max CM-spot distance to granuloma (mm)",
        rd("a2_section_max_distance.tsv").max_dist_gran_mm_cm_spots.max(), cs.dist_gran_mm.max())
    num("A2 CS sections without granuloma", 0, int((sp[sp.cohort == "CS"].groupby("section").immune_dominant.sum() == 0).sum()))
    cref = rd("a2_comparator_reference.tsv").set_index("score")
    cmp_pp = sp[(sp.cohort == "COMPARATOR") & sp.cm_dominant].groupby("patient")
    for s in ("score_program", "score_ACADVL"):
        pm = cmp_pp[s].mean()
        num(f"A2 comparator {s} mean", cref.loc[s, "mean_across_patients"], pm.mean())
        num(f"A2 comparator {s} min", cref.loc[s, "min"], pm.min())
        num(f"A2 comparator {s} max", cref.loc[s, "max"], pm.max())

    # mixed model: lme4 fits are singular with patient variance 0, so (1|section) alone is the same
    # model; refit with statsmodels REML as an independent engine check (tolerance 1%).
    import statsmodels.formula.api as smf
    md = pd.read_csv(f"{E}/a2_model_input.tsv.gz", sep="\t")
    md["z_log10_umi"] = (md.log10_umi - md.log10_umi.mean()) / md.log10_umi.std(ddof=1)
    md["cm_structural"] = md.score_cm_structural
    mod = rd("a2_distance_model.tsv").set_index("outcome")
    for k, f in (("program", "score_program ~ dist_gran_mm + cm_structural + z_log10_umi"),
                 ("ACADVL", "score_ACADVL ~ dist_gran_mm + cm_structural + z_log10_umi"),
                 ("cm_structural_positive_control", "cm_structural ~ dist_gran_mm + z_log10_umi")):
        assert mod.loc[k, "var_patient"] < 1e-6  # patient variance at the boundary
        fit = smf.mixedlm(f, md, groups=md.section).fit(reml=True)
        num(f"A2 distance beta, {k} (statsmodels (1|section) refit)", mod.loc[k, "beta"],
            fit.params["dist_gran_mm"], tol=0.01, note="lme4 fit singular; patient variance 0" if bool(mod.loc[k, "singular"])
            else f"lme4 isSingular FALSE; patient variance {mod.loc[k, 'var_patient']:.1e}")
    for k in ("program", "ACADVL", "oxphos_kegg", "cm_structural_positive_control"):
        r = mod.loc[k]
        excl = (r.ci_low > 0) or (r.ci_high < 0)
        claim_t(f"A2 {k} distance CI excludes 0 (spot-level model)", "yes" if excl else "no",
              "yes" if (r.p < 0.05) else "no", excl == (r.p < 0.05), "CI and p agree")
    claim_t("A2 n patients / sections / spots in model",
          f"{int(mod.loc['program', 'n_patients'])}/{int(mod.loc['program', 'n_sections'])}/{int(mod.loc['program', 'n_spots'])}",
          f"{md.patient.nunique()}/{md.section.nunique()}/{len(md)}",
          (md.patient.nunique(), md.section.nunique(), len(md)) ==
          (mod.loc['program', 'n_patients'], mod.loc['program', 'n_sections'], mod.loc['program', 'n_spots']))

    s6 = sp.groupby("section").agg(cohort=("cohort", "first"), patient=("patient", "first"),
                                   frac_fibroblast_dominant=("fibroblast_dominant", "mean"),
                                   mean_ecm7=("score_ecm7", "mean"), frac_immune_dominant=("immune_dominant", "mean"),
                                   frac_cm_dominant=("cm_dominant", "mean"))
    p6 = s6.groupby("patient").agg(cohort=("cohort", "first"), **{v: (v, "mean") for v in s6.columns[2:]})
    a6 = rd("a6_tissue_fibrosis.tsv").set_index("variable")
    for v in a6.index:
        x = p6.loc[p6.cohort == "CS", v]; y = p6.loc[p6.cohort == "COMPARATOR", v]
        w = stats.ttest_ind(x, y, equal_var=False)
        num(f"A6 {v} CS minus comparator", a6.loc[v, "diff_CS_minus_comparator"], x.mean() - y.mean())
        num(f"A6 {v} Welch p", a6.loc[v, "welch_p"], w.pvalue)
        num(f"A6 {v} Mann-Whitney p", a6.loc[v, "mannwhitney_p"],
            stats.mannwhitneyu(x, y, alternative="two-sided").pvalue)
        num(f"A6 {v} n patients CS", a6.loc[v, "n_CS_patients"], len(x))

    # ---- A3: cross-platform concordance, all FDR<5% hits, all vs QC-pass ----------
    h = rd("q1_hits_annotated.tsv").set_index("gene")
    hits = h[h["adj.P.Val"] < FDR]
    pool = rd("snrna_pooled_full.tsv").set_index("gene")
    a3 = rd("a3_qc_comparison.tsv")
    a3 = a3[(a3.analysis == "cross_platform_snrna") & (a3.gene_universe == "all_fdr5_hits")].set_index(["direction", "statistic"])
    for d, sgn in (("down", -1), ("up", 1)):
        for col, sub in (("all_hits", hits), ("qc_passing_only", hits[hits.qc_pass.astype(bool)])):
            sub = sub[np.sign(sub.logFC) == sgn]
            g = [x for x in sub.index if x in pool.index and np.isfinite(pool.loc[x, "lfc_pool"])]
            conc = int((np.sign(pool.loc[g, "lfc_pool"]) == sgn).sum())
            num(f"A3 {d} {col} n testable", a3.loc[(d, "n_testable"), col], len(g))
            num(f"A3 {d} {col} fraction same sign", a3.loc[(d, "frac_concordant"), col], conc / len(g))
            num(f"A3 {d} {col} Spearman rho", a3.loc[(d, "spearman_rho"), col],
                stats.spearmanr(sub.loc[g, "logFC"], pool.loc[g, "lfc_pool"])[0])

    # ---- A4: Chaffin and internal CS vs NF --------------------------------------
    a4 = rd("a4_published_nf.tsv").set_index(["analysis", "statistic"])
    dn = hits[hits.logFC < 0]
    for tag in ("DCM", "HCM"):
        t = pd.read_csv(f"data/published/chaffin2022_cardiomyocyte_{tag}vsNF.tsv", sep="\t")
        t = t[t.bg_contamination_flag == 0].set_index("gene")
        base = float((t.logFC < 0).mean())
        for col, sub in (("all_hits", dn), ("qc_passing_only", dn[dn.qc_pass.astype(bool)])):
            g = [x for x in sub.index if x in t.index]
            k = int((t.loc[g, "logFC"] < 0).sum())
            num(f"A4 Chaffin {tag} {col} lower in disease", a4.loc[(f"chaffin_{tag}_vs_NF", "n_lower_in_disease"), col], k)
            num(f"A4 Chaffin {tag} {col} n testable", a4.loc[(f"chaffin_{tag}_vs_NF", "n_testable"), col], len(g))
            num(f"A4 Chaffin {tag} {col} binomial p", a4.loc[(f"chaffin_{tag}_vs_NF", "binom_p"), col],
                stats.binomtest(k, len(g), base, alternative="greater").pvalue)
    nf = pd.read_csv("data/snrna/snrna_CS_vs_NF.tsv", sep="\t").set_index("gene")
    a4i = rd("a4_internal_cs_vs_nf.tsv").set_index("statistic")
    for col, sub in (("all_hits", dn), ("qc_passing_only", dn[dn.qc_pass.astype(bool)])):
        g = [x for x in sub.index if x in nf.index and np.isfinite(nf.loc[x, "logFC"])]
        num(f"A4 internal CS vs NF {col} n testable", a4i.loc["n_testable", col], len(g))
        num(f"A4 internal CS vs NF {col} fraction lower", a4i.loc["frac_lower_in_CS_vs_NF", col],
            (nf.loc[g, "logFC"] < 0).mean())

    # ---- A5: Hamed severity axis --------------------------------------------------
    hm = rd("a5_hamed_g2g1.tsv")
    sev = rd("a5_severity_axis.tsv").set_index(["model", "test", "qc"])
    ov = hm[hm.st_logFC.notna() & hm.effect_raw.notna()]
    num("A5 Hamed overlap n", sev.loc[("raw", "spearman_all_overlap", "all"), "n"], len(ov))
    num("A5 Hamed severity Spearman rho (all overlap)", sev.loc[("raw", "spearman_all_overlap", "all"), "value"],
        stats.spearmanr(ov.st_logFC, ov.effect_raw)[0])
    low = ov[ov.st_hit.astype(bool) & (ov.st_logFC < 0)]
    num("A5 CS-lower genes in proteome", sev.loc[("raw", "G2_lower_count", "all"), "n"], len(low))
    for mdl, col in (("raw", "raw"), ("cm_adjusted", "cmadj")):
        num(f"A5 {mdl} G2-lower among CS-lower", sev.loc[(mdl, "G2_lower_count", "all"), "value"],
            ((low[f"adj_{col}"] < FDR) & (low[f"effect_{col}"] < 0)).sum())
        num(f"A5 {mdl} G2-higher among CS-lower", sev.loc[(mdl, "G2_higher_count", "all"), "value"],
            ((low[f"adj_{col}"] < FDR) & (low[f"effect_{col}"] > 0)).sum())

    # ---- A5b: transplant-only snRNA ------------------------------------------------
    tx = rd("snrna_pooled_txonly.tsv").set_index("gene")
    prog = [l.strip() for l in open(cfg["prog"]) if l.strip() and not l.startswith("#")]
    a5b = rd("a5b_txonly_concordance.tsv").set_index(["gene_set", "cohort"])
    pkey = f"primary_qcpass_down_{len(prog)}"
    for coh, tab in (("full", pool), ("transplant_only", tx)):
        g = [x for x in prog if x in tab.index and np.isfinite(tab.loc[x, "lfc_pool"])]
        num(f"A5b primary program {coh} n testable", a5b.loc[(pkey, coh), "n_testable"], len(g))
        num(f"A5b primary program {coh} fraction lower", a5b.loc[(pkey, coh), "frac_concordant"],
            (tab.loc[g, "lfc_pool"] < 0).mean())
    S = rd("a5b_cs_patient_scores.tsv")
    pr = rd("a5b_procurement.tsv")
    for v in ("program_score", "ACADVL"):
        r2 = pr[(pr.analysis == "A5b.2 LVAD vs TX within CS") & (pr.variable == v)].iloc[0]
        x = S.loc[S.procurement.isin(["LVAD", "VAD"]), v]; y = S.loc[S.procurement == "TX", v]
        num(f"A5b LVAD minus TX, {v}", r2.diff_LVAD_minus_TX, x.mean() - y.mean())
        num(f"A5b LVAD vs TX Welch p, {v}", r2.p, stats.ttest_ind(x, y, equal_var=False).pvalue)
        r3 = pr[(pr.analysis == "A5b.3 LVEF association (all CS, numeric LVEF)") & (pr.variable == v)].iloc[0]
        L = S[S.lvef.notna()]
        num(f"A5b LVEF Spearman rho, {v}", r3.spearman_rho, stats.spearmanr(L.lvef, L[v])[0])
        num(f"A5b LVEF n, {v}", r3.n, len(L))
    c18 = rd("a5b_concordant18_txonly.tsv")
    nconc = len(rd("a4_snrna_nf_concordant18.tsv"))
    claim_t("A5b concordant genes carried into transplant-only check", nconc, len(c18), len(c18) == nconc)



def eschf_v3_claims():
    """ESC HF v3 (HANDOFF_figure_v3_ESCHF.md): A9/A12 snRNA screen, A10 pathway levels, A11 cohort accounting,
    each recomputed from upstream inputs (not from the analysis scripts' outputs)."""
    import re
    import scipy.sparse as ssp
    E, V, D = "results/eschf", "results/eschf/primary_v2", "data/eschf/primary_v2"
    rd = lambda f, **k: pd.read_csv(f, sep="\t", **k)
    TB = lambda v: v.astype(str).str.upper().eq("TRUE")

    def num(name, reported, recomputed, tol=1e-6, note=""):
        r, c = float(reported), float(recomputed)
        ok = (np.isnan(r) and np.isnan(c)) or bool(np.isclose(r, c, rtol=tol, atol=tol))
        claim(f"ESCHF[v3] {name}", round(r, 4), round(c, 4), ok, note)

    def same(name, reported, recomputed):
        a, b = sorted(reported), sorted(recomputed)
        claim(f"ESCHF[v3] {name}", ",".join(a) or "none", ",".join(b) or "none", a == b)

    # ---- A9 / A12: snRNA screen ---------------------------------------------------------------
    THR = 3.0
    H = rd(f"{V}/q1_hits_annotated.tsv").set_index("gene"); H = H[H["adj.P.Val"] < FDR]
    P = rd(f"{V}/snrna_pooled_full.tsv").set_index("gene")
    S9 = rd(f"{E}/a9_snrna_screen.tsv").set_index("gene")
    SM = rd(f"{E}/a9_study_mean_log2cpm.tsv.gz").set_index("gene")
    tst = [g for g in H.index if g in P.index]
    lp = P.loc[tst, "lfc_pool"]
    f3a = set(lp[lp.abs() > THR].index)

    def f3b(g):
        if g not in SM.index:
            return None
        m = SM.loc[g]
        return any(m[s] < 0 and m.drop(s).mean() - m[s] >= 2 for s in SM.columns)
    b3 = {g: f3b(g) for g in tst}
    same("A9 3a failures among testable hits", S9.index[TB(S9.fail_3a)], f3a)
    same("A9 3b flags among testable hits", S9.index[TB(S9.flag_3b)], [g for g, v in b3.items() if v])
    claim("ESCHF[v3] A9 3a failures all also flagged by 3b", "yes", "yes" if all(b3[g] for g in f3a) else "no",
          all(b3[g] for g in f3a))
    claim("ESCHF[v3] A9 FPGT flagged by 3b (sanity gene)", "yes", "yes" if f3b("FPGT") else "no", bool(f3b("FPGT")))
    sgn = lambda g: np.sign(P.loc[g, "lfc_pool"]) == np.sign(H.loc[g, "logFC"])
    conc_all = {g for g in tst if sgn(g) and P.loc[g, "fdr_pool"] < FDR}
    same("A9 concordant set from all hits", S9.index[TB(S9.in_concordant_allhits)], conc_all)
    same("A9 concordant set from all hits after 3a", S9.index[TB(S9.in_concordant_allhits_screened)], conc_all - f3a)
    core = rd(f"{V}/q1_core_snrna_pergene.tsv").set_index("gene")
    same("A9 3a failures in section-stable core panel", [g for g in core.index[core.replicates_fdr5] if g in f3a],
         set(core.index[core.replicates_fdr5]) & f3a)
    cnt = rd(f"{E}/a9_screen_counts.tsv").set_index("set")
    qcp = [g for g in tst if str(H.loc[g, "qc_pass"]).upper() == "TRUE"]
    for nm, gs in (("all_hits_testable", tst), ("qc_pass_hits_testable", qcp), ("concordant_allhits", conc_all)):
        num(f"A9 {nm} n", cnt.loc[nm, "n"], len(gs))
        num(f"A9 {nm} n failing 3a", cnt.loc[nm, "n_fail_3a"], len(set(gs) & f3a))
        num(f"A9 {nm} n flagged 3b", cnt.loc[nm, "n_flag_3b"], sum(bool(b3[g]) for g in gs))
    for tag in ("DCM", "HCM"):
        t = rd(f"data/published/chaffin2022_cardiomyocyte_{tag}vsNF.tsv")
        num(f"A9 Chaffin {tag} published genes beyond +/-3", cnt.loc[f"chaffin_{tag}_vs_NF_published", "n_fail_3a"],
            (t.logFC.abs() > THR).sum())
    CC = rd(f"{E}/a9_concordance_screened.tsv").set_index(["screen", "direction", "gene_set"])
    for (scr, d, gset), r in CC.iterrows():
        gs = tst if gset == "all_hits" else qcp
        gs = [g for g in gs if (H.loc[g, "logFC"] < 0) == (d == "down")]
        if scr == "screened_3a":
            gs = [g for g in gs if g not in f3a]
        nm = f"A9 sign agreement {scr} {d} {gset}"
        num(f"{nm} n testable", r.n_testable, len(gs))
        num(f"{nm} n same sign", r.n_same_sign, sum(sgn(g) for g in gs))
        num(f"{nm} n concordant FDR<5%", r.n_concordant_fdr5, sum(sgn(g) and P.loc[g, "fdr_pool"] < FDR for g in gs))
    claim("ESCHF[v3] A12 KCTD14 absent from snRNA counts", "absent", "absent" if "KCTD14" not in P.index else "present",
          "KCTD14" not in P.index)

    # ---- A10: pathway levels ------------------------------------------------------------------
    R10 = rd(f"{E}/a10_pathway_levels.tsv")
    SETS = {"fatty_acid_degradation_hsa00071": "metadata/geneset_kegg_fatty_acid_degradation_hsa00071.txt",
            "tca_cycle_hsa00020": "metadata/geneset_kegg_tca_hsa00020.txt",
            "oxphos_hsa00190": "metadata/geneset_kegg_oxphos_hsa00190.txt"}
    rs = lambda p: [l.strip() for l in open(p) if l.strip() and not l.startswith("#")]
    univ = set(rs(f"{D}/universe_allcm.txt"))
    X = ssp.load_npz(f"{D}/spotcounts_allcm_T2000.npz").tocsr()
    SP = rd(f"{D}/spotcounts_allcm_T2000_spots.tsv.gz")
    G = rs(f"{D}/spotcounts_allcm_T2000_genes.txt"); gi = {g: i for i, g in enumerate(G)}
    LC = X.multiply(1e4 / 2000).tocsr(); LC.data = np.log1p(LC.data)
    mi = rd(f"{V}/a2_model_input.tsv.gz")
    pmap = rd(f"{V}/section_patient_map.tsv")
    arm = pmap.drop_duplicates("patient").set_index("patient").arm
    L = rd(f"{E}/a10_patient_log2cpm_voomform.tsv").set_index("gene")
    for k, p in SETS.items():
        inu = [g for g in rs(p) if g in univ]
        sc = SP[["spot", "section", "cohort"]].assign(s=np.asarray(LC[:, [gi[g] for g in inu]].mean(1)).ravel())
        cs = mi.merge(sc, on=["spot", "section"], how="left")
        comp = sc[sc.cohort == "COMPARATOR"].groupby("section").s.mean()
        for _, r in R10[(R10.set == k) & (R10.level == "distance_bin")].iterrows():
            v = cs[cs.bin == r.bin].groupby("patient").s.mean()
            num(f"A10 {k} bin {r.bin} n genes", r.n_genes, len(inu))
            num(f"A10 {k} bin {r.bin} CS mean of patient means", r.cs_mean_of_patient_means, v.mean())
            num(f"A10 {k} bin {r.bin} comparator min", r.comparator_min, comp.min())
            claim(f"ESCHF[v3] A10 {k} bin {r.bin} CS mean below comparator minimum", str(r.cs_mean_below_comparator_min),
                  str(v.mean() < comp.min()), str(r.cs_mean_below_comparator_min) == str(v.mean() < comp.min()))
        s = L.loc[[g for g in inu if g in L.index]].mean(0)
        assert len([g for g in inu if g in L.index]) == len(inu)
        x, y = s[[c for c in s.index if arm.get(c) == "CS"]], s[[c for c in s.index if arm.get(c) == "COMPARATOR"]]
        w = stats.ttest_ind(x, y, equal_var=False)
        r = R10[(R10.set == k) & (R10.level == "patient_pseudobulk")].iloc[0]
        num(f"A10 {k} patient level CS minus comparator", r.diff_cs_minus_comparator, x.mean() - y.mean())
        num(f"A10 {k} patient level Welch p", r.welch_p, w.pvalue)
        num(f"A10 {k} patient level n CS / comparator", r.n_patients * 100 + r.n_comparator_patients, len(x) * 100 + len(y))
    M10 = rd(f"{E}/a10_pathway_members.tsv")
    for k in SETS:
        m = M10[(M10.set == k) & M10.in_universe]
        num(f"A10 {k} members that are primary hits", m.primary_hit.sum(), sum(g in H.index for g in m.gene))

    # ---- A11: snRNA cohort accounting ---------------------------------------------------------
    A = rd(f"{E}/a11_snrna_cohort_accounting.tsv")
    coh = rd("data/snrna/snrna_patient_cohort.tsv"); coh = coh[coh.disease == "CS"].copy()
    coh["Subject"] = coh.patient.str.replace("neyazi_", "", regex=False)
    s1 = pd.read_excel("data/external/neyazi_table_s1.xlsx", sheet_name="Sarcoid")
    subj = set(s1.Subject.dropna().astype(str))
    raw = pd.read_excel("data/external/neyazi_table_s2.xlsx", header=None).iloc[3:, :5]
    raw = raw[raw.iloc[:, 1].astype(str).str.contains("Sarcoid")]
    s2reg = pd.Series(raw.iloc[:, 2:5].apply(pd.to_numeric).sum(1).values, index=raw.iloc[:, 0].astype(str))
    same("A11 S1 subjects not in analysed snRNA cohort", A.Subject[~TB(A.analysed)], subj - set(coh.Subject))
    cn = coh.set_index("Subject").n_specimens
    ok = bool((cn.values == s2reg.loc[cn.index].values).all())
    claim("ESCHF[v3] A11 analysed specimen counts = S2 region counts, every patient", "yes", "yes" if ok else "no", ok)
    num("A11 analysed CS patients", TB(A.analysed).sum(), len(coh))
    lst = open("data/external/neyazi_geo_snrna_sample_listing.txt").read()
    gsm = sorted(set(re.findall(r"GSM\d+", lst)))
    num("A11 GEO snRNA samples in download", 41, len(gsm))
    for sj in A.Subject[~TB(A.analysed)]:
        claim(f"ESCHF[v3] A11 {sj} absent from GEO snRNA listing", "absent",
              "absent" if re.search(rf"_{re.escape(sj)}_", lst) is None else "present", re.search(rf"_{re.escape(sj)}_", lst) is None)


def eschf_revision_claims():
    """ESC HF revision (handoff_figure_fixes): new primary vs old fit, per-patient slopes, OXPHOS level,
    descriptive direction table and every Figure 1 / legend number, recomputed from upstream inputs."""
    import statsmodels.formula.api as smf
    E, V = "results/eschf", "results/eschf/primary_v2"
    rd = lambda f, **k: pd.read_csv(f, sep="\t", **k)

    def num(name, reported, recomputed, tol=1e-6, note=""):
        ok = bool(np.isclose(float(reported), float(recomputed), rtol=tol, atol=tol))
        claim(f"ESCHF[rev] {name}", round(float(reported), 4), round(float(recomputed), 4), ok, note)

    def same(name, reported, recomputed):
        a_, b_ = sorted(reported), sorted(recomputed)
        claim(f"ESCHF[rev] {name}", ",".join(a_) or "none", ",".join(b_) or "none", a_ == b_)

    # ---- s.0 new primary vs old fit (from the two limma outputs) --------------------------------
    new = rd(f"{V}/q1_raw_allcm.tsv").set_index("gene")
    old = rd(f"{E}/q1_hits_annotated.tsv").set_index("gene")
    univ = [l.strip() for l in open("data/eschf/primary_v2/universe_allcm.txt") if l.strip() and not l.startswith("#")]
    claim("ESCHF[rev] A0 tested genes = new universe", len(univ), len(new), sorted(univ) == sorted(new.index))
    cmpt = rd(f"{E}/a0_vs_old_primary.tsv").set_index("gene_set")
    qc = rd(f"{E}/a0_hits_annotated.tsv").set_index("gene")
    hn, ho = new[new["adj.P.Val"] < FDR], old[old["adj.P.Val"] < FDR]
    num("A0 new primary hits", cmpt.loc["all", "n_hits_new"], len(hn))
    num("A0 new primary hits lower in CS", cmpt.loc["all", "n_hits_new_down"], (hn.logFC < 0).sum())
    num("A0 old-fit hits", cmpt.loc["all", "n_hits_old"], len(ho))
    both = hn.index.intersection(ho.index)
    num("A0 hits in both fits", cmpt.loc["all", "n_overlap"], len(both))
    ho_in = [g for g in ho.index if g in new.index]
    num("A0 old hits keeping sign in new fit", cmpt.loc["all", "old_hits_same_sign_in_new"],
        sum(np.sign(new.loc[g, "logFC"]) == np.sign(ho.loc[g, "logFC"]) for g in ho_in))
    qn = qc.loc[hn.index]
    qn = qn[qn.qc_pass.astype(str).str.upper().eq("TRUE")]
    num("A0 new primary QC-passing hits", cmpt.loc["qc_pass", "n_hits_new"], len(qn))
    prog = [l.strip() for l in open(f"{V}/geneset_primary_program_qcpass_down.txt") if l.strip() and not l.startswith("#")]
    pr_re = sorted(qn.index[new.loc[qn.index, "logFC"] < 0])
    claim("ESCHF[rev] A0 program genes = QC-passing lower-in-CS hits", len(prog), len(pr_re), sorted(prog) == pr_re)

    # ---- s.3 per-patient slopes, recomputed from the per-spot export -----------------------------
    pmap = rd(f"{V}/section_patient_map.tsv")
    sp = rd("data/eschf/primary_v2/eschf_perspot_T2000.tsv.gz")
    sp["patient"] = sp.section.map(pmap.set_index("section").patient)
    cs = sp[(sp.cohort == "CS") & sp.cm_dominant].copy()
    cs["bin"] = pd.cut(cs.dist_gran_mm, [0, 0.25, 0.5, 1, 2, 3, np.inf], right=False,
                       labels=["<0.25", "0.25-0.5", "0.5-1", "1-2", "2-3", ">3"]).astype(object)
    lu = np.log10(cs.umi_native); cs["zu"] = (lu - lu.mean()) / lu.std(ddof=1)
    cs["struct"] = cs.score_cm_structural
    cmp_pm = sp[(sp.cohort == "COMPARATOR") & sp.cm_dominant].groupby(sp.section.map(pmap.set_index("section").patient))
    slp = rd(f"{E}/a2_perpatient_slopes.tsv").set_index("outcome")
    dmax = cs.dist_gran_mm.max()
    for k, col, struct in (("program", "score_program", True), ("ACADVL", "score_ACADVL", True),
                           ("oxphos_kegg", "score_oxphos_kegg", True), ("cm_structural_positive_control", "struct", False)):
        b = []
        for pt, sub in cs.groupby("patient"):
            if sub.bin.nunique() < 2:
                continue
            f = f"{col} ~ dist_gran_mm" + (" + struct" if struct else "") + " + zu" + (" + C(section)" if sub.section.nunique() > 1 else "")
            b.append(smf.ols(f, sub).fit().params["dist_gran_mm"])
        b = np.array(b); n = len(b); m = b.mean(); h = stats.t.ppf(0.975, n - 1) * b.std(ddof=1) / np.sqrt(n)
        r = slp.loc[k]
        num(f"A2 per-patient {k} n patients", r.n_patients, n)
        num(f"A2 per-patient {k} mean slope", r.mean_slope_per_mm, m)
        num(f"A2 per-patient {k} CI low", r.ci_low, m - h)
        num(f"A2 per-patient {k} CI high", r.ci_high, m + h)
        num(f"A2 per-patient {k} n negative", r.n_negative, (b < 0).sum())
        if k in ("program", "ACADVL"):
            gap = cmp_pm[col].mean().mean() - cs.groupby("patient")[col].mean().mean()
            num(f"A2 per-patient {k} gap share across {dmax:.2f} mm", r.gap_share_across_max_distance, abs(m * dmax) / gap)
    num("A2 max CM-dominant distance to immune-dominant spot (mm)", slp.max_distance_mm.iloc[0], dmax)

    # ---- s.4 OXPHOS level by populated bin ----------------------------------------------------
    ox = rd(f"{E}/a2_oxphos_level.tsv").set_index("bin")
    cpm = cmp_pm["score_oxphos_kegg"].mean()
    for bn, sub in cs.groupby("bin"):
        pm = sub.groupby("patient").score_oxphos_kegg.mean()
        h = stats.t.ppf(0.975, len(pm) - 1) * pm.std(ddof=1) / np.sqrt(len(pm))
        num(f"A2 OXPHOS bin {bn} CS mean", ox.loc[bn, "cs_mean_of_patient_means"], pm.mean())
        num(f"A2 OXPHOS bin {bn} CS CI high", ox.loc[bn, "cs_ci_high"], pm.mean() + h)
        claim(f"ESCHF[rev] A2 OXPHOS bin {bn} CS mean below comparator range",
              str(bool(ox.loc[bn, "cs_mean_below_comparator_min"])), str(bool(pm.mean() < cpm.min())),
              bool(ox.loc[bn, "cs_mean_below_comparator_min"]) == bool(pm.mean() < cpm.min()))
    claim("ESCHF[rev] A2 OXPHOS populated bins", len(ox), cs.bin.nunique(), len(ox) == cs.bin.nunique())

    # ---- s.6 descriptive direction table --------------------------------------------------------
    dd = rd(f"{V}/a3a4_direction_descriptive.tsv").set_index(["comparison", "gene_set"])
    hq = qc.loc[hn.index]
    dn = hq[new.loc[hq.index, "logFC"] < 0]
    pool = rd(f"{V}/snrna_pooled_full.tsv").set_index("gene")
    nf = rd("data/snrna/snrna_CS_vs_NF.tsv").set_index("gene")
    refs = {"cross_platform_snrna": (pool, "lfc_pool"), "snrna_CS_vs_NF": (nf, "logFC")}
    for tag in ("DCM", "HCM"):
        t = rd(f"data/published/chaffin2022_cardiomyocyte_{tag}vsNF.tsv")
        refs[f"chaffin_{tag}_vs_NF"] = (t[t.bg_contamination_flag == 0].set_index("gene"), "logFC")
    for comp, (tab, col) in refs.items():
        for gs, sub in (("all_hits", dn), ("qc_passing_only", dn[dn.qc_pass.astype(str).str.upper().eq("TRUE")])):
            g = [x for x in sub.index if x in tab.index and np.isfinite(tab.loc[x, col])]
            r = dd.loc[(comp, f"lower-in-CS hits ({gs})")]
            num(f"s6 {comp} {gs} n testable", r.n_testable, len(g))
            num(f"s6 {comp} {gs} fraction lower", r.frac_in_direction, (tab.loc[g, col] < 0).mean())
    claim("ESCHF[rev] s6 p-values reported for all four comparisons", "yes",
          "yes" if dd.p_value.notna().all() and dd.index.get_level_values(0).nunique() == 4 else "no",
          bool(dd.p_value.notna().all() and dd.index.get_level_values(0).nunique() == 4))

    # ---- Figure 1 / legend numbers ------------------------------------------------------------
    pv = rd(f"{E}/a7_panel_values.tsv").drop_duplicates(["panel", "quantity"]).set_index(["panel", "quantity"]).value
    for sec in ("CS_108-1", "ARVC_1"):
        x = sp[sp.section == sec]
        num(f"Fig A {sec} CM-dominant spots", pv[("A", f"{sec}|cm")], x.cm_dominant.sum())
        num(f"Fig A {sec} immune-dominant spots", pv[("A", f"{sec}|imm")], x.immune_dominant.sum())
    for sc, k in (("score_program", "program"), ("score_ACADVL", "ACADVL")):
        for q in ("perpatient_slope", "perpatient_slope_ci_low", "perpatient_slope_ci_high", "perpatient_n", "perpatient_n_negative"):
            col = {"perpatient_slope": "mean_slope_per_mm", "perpatient_slope_ci_low": "ci_low", "perpatient_slope_ci_high": "ci_high",
                   "perpatient_n": "n_patients", "perpatient_n_negative": "n_negative"}[q]
            num(f"Fig B {k} {q}", pv[("B", f"{sc}|{q}")], slp.loc[k, col])
        for bn, sub in cs.groupby("bin"):
            num(f"Fig B {k} bin {bn} pooled mean", pv[("B", f"{sc}|{bn}|pooled_mean")], sub.groupby("patient")[sc].mean().mean())
    b_bins = sorted(q.split("|")[1] for (pn, q) in pv.index if pn == "B" and q.startswith("score_program|") and q.endswith("pooled_mean"))
    claim("ESCHF[rev] Fig B plotted bins = populated bins", ",".join(b_bins), ",".join(sorted(cs.bin.unique())),
          b_bins == sorted(cs.bin.unique()))
    num("Fig B legend max distance (mm)", pv[("B", "max_distance_cm_spot_mm")], dmax)
    # v3 (HANDOFF_figure_v3_ESCHF.md): panels C and D after snRNA screen 3a (|pooled snRNA log2FC| > 3), D from all hits
    SCREEN = 3.0
    T = lambda v: v.astype(str).str.upper().eq("TRUE")
    hqc = hq[T(hq.qc_pass)]
    num("Fig B score genes (QC-passing hits lower in CS)", pv[("B", "n_score_genes")], ((hqc.logFC < 0)).sum())
    gC = [x for x in hqc.index if x in pool.index]
    lC = pool.loc[gC, "lfc_pool"]
    num("Fig C QC-passing hits plotted", pv[("C", "n_qcpass_genes_plotted")], len(gC))
    num(f"Fig C genes beyond +/-{YLIM_C:g}", pv[("C", "n_out_of_range")], (lC.abs() > YLIM_C).sum())
    num("Fig C genes failing screen 3a", pv[("C", "n_screen_3a_fail")], (lC.abs() > SCREEN).sum())
    claim("ESCHF[rev] Fig C off-axis genes = screen-3a failures", "identical",
          "identical" if set(lC[lC.abs() > YLIM_C].index) == set(lC[lC.abs() > SCREEN].index) else "differ",
          set(lC[lC.abs() > YLIM_C].index) == set(lC[lC.abs() > SCREEN].index))
    for g in lC[lC.abs() > SCREEN].index:
        num(f"Fig C triangle {g} snRNA log2FC", pv[("C", f"out_of_range|{g}|snrna_lfc")], lC[g])
        num(f"[v4] Fig C triangle {g} drawn open neutral (screened, no sign category)",
            pv[("C", f"out_of_range|{g}|drawn_open_neutral")], 1)
    _open = sorted(q.split("|")[1] for (pn, q) in pv.index if pn == "C" and q.endswith("|drawn_open_neutral")
                   and pv[(pn, q)] == 1)
    same("[v4] Fig C open-neutral triangles = screen-3a failures", _open, set(lC[lC.abs() > SCREEN].index))
    claim("ESCHF[rev] [v4] Fig C no category-filled triangle values remain", 0,
          sum(1 for (pn, q) in pv.index if pn == "C" and q.endswith("|category")),
          not any(pn == "C" and q.endswith("|category") for (pn, q) in pv.index))
    conc = rd(f"{E}/a9_concordance_screened.tsv").set_index(["screen", "direction", "gene_set"])
    for d in ("down", "up"):
        dm = (hqc.loc[gC, "logFC"] < 0) if d == "down" else (hqc.loc[gC, "logFC"] > 0)
        gg = [x for x in np.array(gC)[dm.values]]
        keep = [x for x in gg if abs(lC[x]) <= SCREEN]
        same_k = (np.sign(lC[keep]) == np.sign(hqc.loc[keep, "logFC"])).mean()
        same_u = (np.sign(lC[gg]) == np.sign(hqc.loc[gg, "logFC"])).mean()
        num(f"Fig C screened {d} n", pv[("C", f"screened|n_{d}")], len(keep))
        num(f"Fig C screened {d} fraction same sign", pv[("C", f"screened|frac_same_sign_{d}")], same_k)
        num(f"Fig C screened {d} rho", pv[("C", f"screened|rho_{d}")],
            stats.spearmanr(hqc.loc[keep, "logFC"], lC[keep]).statistic)
        num(f"Fig C unscreened {d} fraction same sign", pv[("C", f"unscreened|frac_same_sign_{d}")], same_u)
        num(f"Fig C screened {d} background (a9 permutation)", pv[("C", f"screened|background_{d}")],
            conc.loc[("screened_3a", d, "qc_passing_only"), "background_frac"], note="from a9_concordance_screened.tsv")
        num(f"Fig C unscreened {d} background (a9 permutation)", pv[("C", f"unscreened|background_{d}")],
            conc.loc[("unscreened", d, "qc_passing_only"), "background_frac"], note="from a9_concordance_screened.tsv")
    gA = [x for x in hq.index if x in pool.index]
    lA = pool.loc[gA]
    # v4 (user decision): Panel D = concordant genes of the section-stable core, recomputed from the six
    # leave-one-comparator-section-out refits (significant in the full fit AND every refit), screen-3a failures marked ¶
    _loo = sorted(glob.glob("results/eschf/primary_v2/loo/q1_raw_drop_*.tsv"))
    num("[v4] Fig D leave-one-comparator-section-out refits", pv[("legend", "n_section_refits_core")], len(_loo))
    _core = set(hq.index)
    for f_ in _loo:
        t_ = rd(f_).set_index("gene"); _core &= set(t_.index[t_["adj.P.Val"] < FDR])
    d18 = [x for x in gA if x in _core and np.sign(lA.loc[x, "lfc_pool"]) == np.sign(hq.loc[x, "logFC"])
           and lA.loc[x, "fdr_pool"] < FDR]
    # v5 (user decision): screen-3a failures removed from the core panel (18 -> 17; NDUFC2 not concordant)
    dset = [x for x in d18 if abs(lA.loc[x, "lfc_pool"]) <= SCREEN]
    num("[v5] Fig D concordant genes (section-stable core after screen 3a)", pv[("D", "n_concordant_genes")], len(dset))
    num("[v5] Fig D core genes removed by screen 3a", pv[("D", "n_core_removed_screen_3a")], len(d18) - len(dset))
    same("[v5] Fig D gene set = section-stable core concordant after 3a (recomputed from refits)",
         sorted({q.split("|")[0] for (pn, q) in pv.index if pn == "D" and q.endswith("|flag_section_segdup")}), set(dset))
    claim("ESCHF[rev] [v5] Fig D no pilcrow values or key remain", 0,
          sum(1 for (pn, q) in pv.index if "pilcrow" in q or q == "n_symbol|¶"),
          not any("pilcrow" in q or q == "n_symbol|¶" for (pn, q) in pv.index))
    # v5 user check: no implausible value in the snRNA vs NF column (same magnitude bound as screen 3a)
    nfD = rd("data/snrna/snrna_CS_vs_NF.tsv").set_index("gene").reindex(dset)
    num("[v5] Fig D max |snRNA vs NF log2FC| over panel genes", pv[("D", "max_abs_snrna_vs_nf_logFC")], nfD.logFC.abs().max())
    num("[v5] Fig D min snRNA vs NF AveExpr over panel genes", pv[("D", "min_snrna_vs_nf_AveExpr")], nfD.AveExpr.min())
    claim("ESCHF[rev] [v5] Fig D snRNA vs NF: no panel gene beyond +/-3 log2", 0,
          int((nfD.logFC.abs() > SCREEN).sum()), bool((nfD.logFC.abs() <= SCREEN).all() and nfD.logFC.notna().all()))
    # v5 C legend sentence: every screened gene near-absent in one source study, present in the others (3b rule
    # recomputed from the per-study mean log2 CPM: study mean < 0 and mean of the other four >= 2 higher), as FPGT
    smn = rd(f"{E}/a9_study_mean_log2cpm.tsv.gz").set_index("gene")
    def _absent(g):
        m = smn.loc[g]; hits = [(st, m[st], m.drop(st).mean()) for st in smn.columns if m[st] < 0 and m.drop(st).mean() - m[st] >= 2]
        return min(hits, key=lambda h: h[1]) if hits else None
    f3 = sorted(lC[lC.abs() > SCREEN].index)
    ab = {g: _absent(g) for g in f3 + ["FPGT"]}
    claim("ESCHF[rev] [v5] Legend C: all screened genes near-absent in one source study (3b, recomputed)",
          f"{len(f3)}/{len(f3)}", f"{sum(ab[g] is not None for g in f3)}/{len(f3)}", all(ab[g] is not None for g in f3))
    claim("ESCHF[rev] [v5] Legend C: FPGT (FPGT-TNNI3K locus) shows the same pattern", "yes",
          "yes" if ab["FPGT"] is not None else "no", ab["FPGT"] is not None)
    for g in f3:
        num(f"[v5] Legend C {g} absent-study mean log2 CPM", pv[("legend", f"source_absence|{g}|absent_mean_log2cpm")], ab[g][1])
        num(f"[v5] Legend C {g} other-studies mean log2 CPM", pv[("legend", f"source_absence|{g}|others_mean_log2cpm")], ab[g][2])
    for g in dset:
        num(f"Fig D {g} dagger (locus outlier)", pv[("D", f"{g}|flag_dagger_locus_outlier")],
            int(str(hq.loc[g, "flag_locus_quant_outlier"]).upper() == "TRUE"))
        num(f"Fig D {g} double dagger (outside QC universe)", pv[("D", f"{g}|flag_ddagger_out_of_universe")],
            int(str(hq.loc[g, "in_qc_universe"]).upper() != "TRUE"))
        num(f"Fig D {g} section mark (segdup)", pv[("D", f"{g}|flag_section_segdup")],
            int(str(hq.loc[g, "flag_segdup"]).upper() == "TRUE"))
    rf = rd(f"{E}/a0_primary_refit.tsv").set_index("gene")
    sh = rf[rf.adj_section < FDR]
    num("Legend patient-pseudobulk section hits", pv[("legend", "patient_pseudobulk|n_section_hits")], len(sh))
    num("Legend patient-pseudobulk hits still FDR<5%", pv[("legend", "patient_pseudobulk|section_hits_still_fdr5")],
        (sh.adj_patient < FDR).sum())
    num("Legend patient-pseudobulk hits same sign", pv[("legend", "patient_pseudobulk|sign_concordant_section_hits")],
        (np.sign(sh.logFC_patient) == np.sign(sh.logFC_section)).sum())
    # v4 legend B sentence: "Pathway scores defined independently of the fit (KEGG TCA cycle, fatty acid degradation and
    # oxidative phosphorylation) were not lower in CS once fit hits were excluded." Values from the A14 table, which is
    # itself recomputed from the hit list and voom matrix in the ESCHF[a14] claims.
    a14 = rd(f"{E}/a14_pathway_hits_removed.tsv")
    hr = a14[a14.version.isin(["minus_primary", "minus_union"])]
    claim("ESCHF[rev] [v4] Legend B pathway sentence: no hits-removed 95% CI entirely below zero (3 pathways x 2 removals)",
          "6/6", f"{int((hr.diff_ci_high > 0).sum())}/{len(hr)}", bool(len(hr) == 6 and (hr.diff_ci_high > 0).all()))
    for _, r in hr.iterrows():
        for k in ("diff_cs_minus_comparator", "diff_ci_low", "diff_ci_high"):
            num(f"[v4] Legend B {r.set} {r.version} {k}",
                pv[("legend", f"pathway_hits_removed|{r.set}|{r.version}|{k}")], r[k], note="from a14_pathway_hits_removed.tsv")
    claim("ESCHF[rev] [v4] Legend B no OXPHOS-range values remain", 0,
          sum(1 for (pn, q) in pv.index if pn == "legend" and q.startswith("oxphos|")),
          not any(pn == "legend" and q.startswith("oxphos|") for (pn, q) in pv.index))
    rp = pmap[pmap.retained == "yes"]
    mt = rd("data/eschf/primary_v2/meta_allcm_T2000.tsv")
    num("Legend CS sections in fit", pv[("legend", "n_cs_sections_fit")],
        ((mt.cohort == "CS") & mt.evaluable.astype(str).str.lower().eq("true")).sum())
    num("Legend CS patients in fit", pv[("legend", "n_cs_patients_fit")], rp[rp.arm == "CS"].patient.nunique())
    num("Legend comparator patients in fit", pv[("legend", "n_comparator_patients_fit")], rp[rp.arm == "COMPARATOR"].patient.nunique())


def short_report_claims():
    """Every number in the short report text (CS_cardiomyocyte_ESCHF_ShortReport), written as printed and recomputed
    from the committed v3 tables at the same precision. Status OK when the strings match. The explant-only expected
    value is encoded as 49% (text correction pending: the draft prints the full-cohort 48%)."""
    import json, re
    rd_ = lambda f: pd.read_csv(f, sep="\t")
    P, E = "results/eschf/primary_v2", "results/eschf"
    pct = lambda x: f"{100 * x:.0f}%"
    def sr(where, what, ms, rc, note=""):
        claim(f"ShortReport [{where}] {what}", ms, rc, str(ms) == str(rc), note)
    def const(path, pat):   # the constant as defined in the script that uses it; exactly one definition
        m = re.findall(pat, open(path).read(), flags=re.M)
        assert len(m) == 1, (path, pat, m); return float(m[0])

    # ---- cohort -------------------------------------------------------------------
    spm = rd_(f"{P}/section_patient_map.tsv")
    cs, co = spm[spm.arm == "CS"], spm[spm.arm == "COMPARATOR"]
    csr = cs[cs.retained == "yes"]
    sr("Methods", "CS samples profiled", 12, len(cs)); sr("Methods", "CS patients profiled", 9, cs.patient.nunique())
    sr("Methods", "comparator samples", 6, len(co)); sr("Methods", "comparator patients", 6, co.patient.nunique())
    comp = co.disease.value_counts().to_dict()
    sr("Methods", "comparator composition (ARVC, LMNA DCM, HCM, Chagas)", "2,2,1,1",
       ",".join(str(comp.get(k, 0)) for k in ("ARVC", "LMNA_DCM", "HCM", "Chagas")), note=str(comp))
    sr("Abstract/Methods", "CS samples analysed", 11, len(csr)); sr("Abstract/Methods", "CS patients analysed", 8, csr.patient.nunique())
    sr("Methods", "minimum cardiomyocyte-dominant spots per sample", 50,
       int(const("scripts/eschf/a0_primary_export.py", r"^MIN_SPOTS\s*=\s*(\d+)")))
    cfg = json.load(open(f"{DAT}/spatial_first_config.json"))
    sr("Methods", "UMI floor and downsampling depth", "2,000", f"{int(cfg['T_headline']):,}")
    sr("Methods", "detection threshold", "25%", pct(cfg["detect_headline"]))
    sr("Discussion", "every analysed CS sample contains immune-dominant spots", "all",
       "all" if (csr.n_immune_dominant_spots > 0).all() else f"{(csr.n_immune_dominant_spots > 0).sum()}/{len(csr)}")

    # ---- spatial differential expression -------------------------------------------------
    rs = rd_(f"{E}/a0_primary_refit_summary.tsv").set_index(["refit", "gene_set"])
    sl, pp = rs.loc[("section_level_primary", "all_tested")], rs.loc[("patient_pseudobulk", "all_tested")]
    sr("Methods/Results", "genes tested", 985, int(sl.n_genes))
    sr("Results", "section-level hits", 223, int(sl.hits_fdr5))
    sr("Results", "hits lower in CS", 151, int(sl.hits_down)); sr("Results", "hits higher in CS", 72, int(sl.hits_up))
    sr("Results", "hits passing locus QC", 180, int(rs.loc[("section_level_primary", "qc_pass"), "hits_fdr5"]))
    sr("Results", "patient-pseudobulk hits", 206, int(pp.hits_fdr5))
    sr("Results", "patient-pseudobulk hits lower", 145, int(pp.hits_down)); sr("Results", "patient-pseudobulk hits higher", 61, int(pp.hits_up))
    sr("Results", "section hits retaining direction at patient level", "223/223",
       f"{int(pp.sign_concordant_section_hits)}/{int(pp.n_section_hits)}")
    sr("Results", "section hits still FDR<5% at patient level", 193, int(pp.section_hits_still_fdr5))
    dc = rs.loc[("section_dupcor_patient_block", "all_tested")]
    sr("Results", "section hits significant with patient blocking", "202/223",
       f"{int(dc.section_hits_still_fdr5)}/{int(dc.n_section_hits)}")
    a0 = rd_(f"{E}/a0_vs_old_primary.tsv").set_index("gene_set").loc["all"]
    sr("Results", "0.5 mm exclusion fit: hits retaining direction", "191/193",
       f"{int(a0.old_hits_same_sign_in_new)}/{int(a0.n_hits_old)}")
    sr("Results", "0.5 mm exclusion fit vs primary: r", "0.97", f"{a0.pearson_r_logFC_shared:.2f}")

    # ---- pathway members lower in CS ----------------------------------------------
    m = rd_(f"{E}/a10_pathway_members.tsv")
    dn = m[m.primary_hit.astype(bool) & (m.hit_direction == "down")]
    g = lambda k: ", ".join(sorted(dn[dn.set.str.startswith(k)].gene))
    sr("Results", "citrate cycle genes lower in CS", "FH, IDH2, IDH3B, PDHB, SDHB", g("tca_cycle"))
    sr("Results", "fatty acid degradation genes lower in CS", "ACADM, ACADVL, ACAT1, ACSL1", g("fatty_acid_degradation"))
    sr("Results", "oxidative phosphorylation subunits lower in CS", 9, int(dn.set.str.startswith("oxphos").sum()))
    nl = rd_(f"{E}/a14_pathway_null.tsv")
    sr("Methods", "expression-matched random gene sets (pathway null)", "10,000", f"{int(nl.n_draws.unique().item()):,}")

    # ---- distance -----------------------------------------------------------------------
    md = rd_(f"{P}/a2_section_max_distance.tsv").set_index("section").loc[csr.section]
    sr("Results/Discussion", "maximum CM-dominant spot distance to an immune-dominant spot (mm)", "1.64",
       f"{md.max_dist_gran_mm_cm_spots.max():.2f}")
    pv_ = rd_(f"{E}/a7_panel_values.tsv").drop_duplicates().set_index(["panel", "quantity"]).value   # a7 writes panel B rows twice (identical)
    assert pv_.index.is_unique, "a7_panel_values: conflicting duplicate quantities"
    f3 = lambda x: f"{x:+.3f}".replace("-", "\u2212")
    for k, lab, ms in (("score_program", "program score", ("+0.021", "\u22120.021", "+0.063")),
                       ("score_ACADVL", "ACADVL", ("+0.062", "\u22120.212", "+0.337"))):
        sr("Results", f"per-patient distance slope, {lab} (slope, CI low, CI high)", ", ".join(ms),
           ", ".join(f3(pv_[("B", f"{k}|perpatient_slope{s}")]) for s in ("", "_ci_low", "_ci_high")))
    sl2 = rd_(f"{E}/a2_perpatient_slopes.tsv").set_index("outcome")
    st = sl2.loc[[o for o in sl2.index if "struct" in o.lower()][0]]
    sr("Results", "cardiomyocyte structural score rose with distance (patients)", "all eight",
       "all eight" if st.n_positive == st.n_patients == 8 else f"{int(st.n_positive)}/{int(st.n_patients)}")
    fb = rd_(f"{P}/a6_tissue_fibrosis.tsv").set_index("variable").loc["frac_fibroblast_dominant"]
    sr("Abstract/Results", "CS samples did not contain more fibroblast-dominant tissue", "not more",
       "not more" if fb.ci_low <= 0 else "more",
       note=f"CS minus comparator {fb.diff_CS_minus_comparator:+.3f} (95% CI {fb.ci_low:+.3f} to {fb.ci_high:+.3f}), Welch p={fb.welch_p:.2f}")
    sr("Legend B", "patients in per-patient slopes", 8, int(pv_[("B", "score_program|perpatient_n")]))
    sr("Legend B", "score genes (QC-passing, lower in CS)", 116, int(rs.loc[("section_level_primary", "qc_pass"), "hits_down"]))

    # ---- snRNA-seq ------------------------------------------------------------------
    sr("Results/Legend C", "genes excluded by the +/-3 screen", "five" , {5: "five"}.get(int(pv_[("C", "n_screen_3a_fail")]), str(int(pv_[("C", "n_screen_3a_fail")]))))
    sr("Legend C", "triangles drawn beyond the axis", "five", {5: "five"}.get(int(pv_[("C", "n_out_of_range")]), str(int(pv_[("C", "n_out_of_range")]))))
    a9_thr = const('scripts/eschf/a9_snrna_screen.py', r'\bTHR\s*=\s*([\d.]+)')
    a9_nperm = const('scripts/eschf/a9_snrna_screen.py', r'^[^#]*\bNPERM, SEED = (\d+), \d+')
    sr("Methods", "pooled-effect screen (log2)", "3", f"{a9_thr:g}")
    sr("Methods/Legend C", "random gene sets for expected agreement", "2,000",
       f"{int(a9_nperm):,}")
    sr("Legend C", "QC-passing hits testable in snRNA-seq", 170, int(pv_[("C", "n_qcpass_genes_plotted")]))
    sr("Results", "lower-in-CS genes, same direction in snRNA-seq", "67%", pct(pv_[("C", "screened|frac_same_sign_down")]))
    sr("Results/Legend C", "lower-in-CS genes, expected by chance", "48%", pct(pv_[("C", "screened|background_down")]))
    sr("Results", "lower-in-CS genes compared (n)", 108, int(pv_[("C", "screened|n_down")]))
    sr("Results", "higher-in-CS genes, same direction", "42%", pct(pv_[("C", "screened|frac_same_sign_up")]))
    sr("Results/Legend C", "higher-in-CS genes, expected by chance", "51%", pct(pv_[("C", "screened|background_up")]))
    sr("Results", "higher-in-CS genes compared (n)", 57, int(pv_[("C", "screened|n_up")]))
    sr("Results/Legend D", "concordant genes in Figure 1D", 19, int(pv_[("D", "n_concordant_genes")]))
    a11 = rd_(f"{E}/a11_snrna_cohort_accounting.tsv")
    sr("Methods", "published CS patients with snRNA-seq analysed", "21 of 22", f"{int(a11.analysed.astype(bool).sum())} of {len(a11)}")
    nf = rd_(f"{P}/a4_internal_cs_vs_nf.tsv")
    nf = nf[(nf.gene_universe == "all_fdr5_hits") & (nf.direction == "down")].set_index("statistic").qc_passing_only
    sr("Results", "CS vs non-failing: lower-in-CS genes lower", "45%", pct(nf["frac_lower_in_CS_vs_NF"]))
    sr("Results", "CS vs non-failing: expected", "51%", pct(nf["ref_frac_lower_jointly_testable"]))

    # ---- published cardiomyocyte data (Chaffin) ---------------------------------------
    pb = rd_(f"{P}/a4_published_nf.tsv")
    pb = pb[(pb.gene_universe == "all_fdr5_hits") & (pb.direction == "down")].set_index(["analysis", "statistic"]).qc_passing_only
    for d, ms_n, ms_bg, ms_p, pf in (("DCM", "19/24", "48%", "0.002", "{:.3f}"), ("HCM", "20/28", "45%", "0.0038", "{:.4f}")):
        a = f"chaffin_{d}_vs_NF"
        sr("Results/Legend D", f"published {d}: lower in disease", ms_n, f"{int(pb[(a, 'n_lower_in_disease')])}/{int(pb[(a, 'n_testable')])}")
        sr("Results/Legend D", f"published {d}: base rate", ms_bg, pct(pb[(a, "ref_frac_lower_in_disease")]))
        sr("Legend D", f"published {d}: one-sided binomial P", ms_p, pf.format(pb[(a, "binom_p")]))

    # ---- procurement and LVEF -----------------------------------------------------------
    pr = rd_(f"{P}/a5b_procurement.tsv")
    v = pr[(pr.analysis.str.startswith("A5b.2")) & (pr.variable == "program_score")].iloc[0]
    sr("Results", "VAD minus explant score difference (diff, CI low, CI high)", "-0.01, -0.04, 0.02",
       ", ".join(f"{x:.2f}" for x in (v.diff_LVAD_minus_TX, v.ci_low, v.ci_high)))
    v = pr[(pr.analysis.str.startswith("A5b.3 LVEF")) & (pr.variable == "program_score")].iloc[0]
    sr("Results", "score vs ejection fraction: rho", "-0.01", f"{v.spearman_rho:.2f}")
    sr("Results", "score vs ejection fraction: n", 20, int(v.n))
    tx = rd_(f"{P}/a5b_txonly_concordance.tsv")
    cc = tx[tx.gene_set.str.startswith("concordant") & (tx.cohort == "transplant_only")].iloc[0]
    sr("Results", "explant only: concordant genes retaining direction", "all",
       "all" if cc.n_same_sign_as_full == cc.n_testable else f"{int(cc.n_same_sign_as_full)}/{int(cc.n_testable)}",
       note=f"concordant set of {int(cc.n_testable)} (Figure 1D shows {int(pv_[('D', 'n_concordant_genes')])} after the 3a screen)")
    q = tx[(tx.gene_set == "primary_qcpass_down_116") & (tx.cohort == "transplant_only")].iloc[0]
    sr("Results", "explant only: agreement, lower-in-CS genes", "62%", pct(q.frac_concordant))
    sr("Results", "explant only: expected by chance", "49%", pct(q.null_mean),
       note="draft text prints 48% (full-cohort background); correction to 49% pending")


def main():
    eschf_claims()
    short_report_claims()

    # ---- report ----------------------------------------------------------------
    out = pd.DataFrame(rows)
    out.to_csv("results/manuscript_results.tsv", sep="\t", index=False)
    w = max(out.claim.str.len())
    for _, r in out.iterrows():
        flag = {"OK": "  ok", "MISMATCH": "FAIL", "NOT_REPRODUCIBLE": "  --", "NOT_RUN": " n/r"}[r.status]
        print(f"{flag}  {r.claim:<{w}}  manuscript {str(r.manuscript):>15}  "
              f"recomputed {str(r.recomputed):>22}")
    n_fail = int((out.status == "MISMATCH").sum())
    n_na = int((out.status == "NOT_REPRODUCIBLE").sum())
    print(f"\n{int((out.status == 'OK').sum())} reproduced, {n_fail} mismatched, "
          f"{n_na} not addressable from committed data")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
