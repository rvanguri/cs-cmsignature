#!/usr/bin/env python3
"""Section-stability of the spatial contrast, and its cross-platform comparison.

Three things happen here, in order.

1. Leave-one-out stability. The Q1 contrast is refit with one comparator section
   dropped (six refits) and with one comparator disease dropped (four refits; the two
   single-section diseases reuse the corresponding section refit). Each refit is
   summarised against the full-cohort result: how many genes stay significant, whether
   any effect changes sign, and the rank agreement over the full tested set.
   The *section-stable core* is the set of genes significant in the full cohort AND in
   all six section refits, i.e. no single comparator section carries the result.

2. Single-nucleus pooling. The two independent single-nucleus contrasts (CS vs DCM and
   CS vs ARVC) are combined by inverse-variance meta-analysis, with the standard error
   taken as |logFC / t| from the limma fit. This is the second platform; it shares no
   samples with the spatial cohort.

3. Cross-platform concordance. For the stable-down and stable-up halves of the core,
   sign concordance against the pooled single-nucleus estimate is tested two ways: a
   binomial test against 0.5, and an expression-decile-matched permutation null, which
   is the relevant comparison because sign concordance is not 50/50 at random in a
   cohort where one arm is systematically lower-expressing.

Inputs
  data/spatial_first/meta_T2000.tsv
  data/snrna/snrna_CS_vs_DCM.tsv, snrna_CS_vs_ARVC.tsv
  results/spatial_first/q1_raw_T2000_d25.tsv
  results/spatial_first/loo/q1_raw_drop_*.tsv, q1_raw_dropdis_*.tsv   (step 02)

Outputs (all under results/spatial_first/)
  q1_leave_one_comparator_out.tsv   per-section refit summary
  q1_leave_one_disease_out.tsv      per-disease refit summary
  q1_stable_down_115.tsv            the stable-down core, annotated with snRNA pooling
  q1_core_snrna_pergene.tsv         per-gene cross-platform table (figure panel A input)
  q1_core_snrna_concordance.tsv     concordance tests for the two core halves
"""
import argparse

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

DATA = "data/spatial_first"
SNRNA = "data/snrna"
RES = "results/spatial_first"
FDR = 0.05
NPERM = 2000
SEED = 0

# Comparator diseases contributing more than one section get their own disease-level
# refit; the single-section diseases are covered by the section-level refit.
DISEASE_SECTIONS = {"ARVC": ["ARVC_1", "ARVC_2"], "LMNA_DCM": ["DC_LMNA", "LMNA"]}


def pool_snrna():
    """Inverse-variance pooling of the two single-nucleus comparator contrasts."""
    sn = {}
    for tag in ("DCM", "ARVC"):
        t = pd.read_csv(f"{SNRNA}/snrna_CS_vs_{tag}.tsv", sep="\t").set_index("gene")
        t["se"] = (t.logFC / t.t).abs()
        sn[tag] = t
    P = pd.DataFrame({
        "lfc_DCM": sn["DCM"].logFC, "se_DCM": sn["DCM"].se, "fdr_DCM": sn["DCM"]["adj.P.Val"],
        "lfc_ARVC": sn["ARVC"].logFC, "se_ARVC": sn["ARVC"].se, "fdr_ARVC": sn["ARVC"]["adj.P.Val"],
        "AveExpr": sn["DCM"].AveExpr,
    }).dropna(subset=["lfc_DCM", "lfc_ARVC", "se_DCM", "se_ARVC"])
    wD, wA = 1 / P.se_DCM ** 2, 1 / P.se_ARVC ** 2
    P["lfc_pool"] = (wD * P.lfc_DCM + wA * P.lfc_ARVC) / (wD + wA)
    P["se_pool"] = np.sqrt(1 / (wD + wA))
    P["z_pool"] = P.lfc_pool / P.se_pool
    P["p_pool"] = 2 * stats.norm.sf(P.z_pool.abs())
    P["fdr_pool"] = multipletests(P.p_pool, method="fdr_bh")[1]
    P["expr_dec"] = pd.qcut(P.AveExpr, 10, labels=False, duplicates="drop")
    return P


def gene_length_scan():
    """Per-contrast regression of snRNA-seq log2FC on log10 gene span.

    Documents why CS vs HCM is excluded from the pooled comparator: its effect sizes
    track gene length and most genes are called differential.
    """
    span = pd.read_csv(f"{DATA}/gene_span.tsv", sep="\t", index_col=0)["span_kb"]
    rows = []
    for tag in ("HCM", "DCM", "ARVC", "NF"):
        t = pd.read_csv(f"{SNRNA}/snrna_CS_vs_{tag}.tsv", sep="\t").set_index("gene")
        t = t.join(span, how="inner").dropna(subset=["span_kb"])
        t = t[t.span_kb > 0]
        lr = stats.linregress(np.log10(t.span_kb), t.logFC)
        sig = t["adj.P.Val"] < FDR
        rows.append(dict(
            contrast=f"CS_vs_{tag}", n_genes=len(t),
            slope_logFC_per_log10span=round(lr.slope, 3), p_slope=lr.pvalue,
            r=round(lr.rvalue, 3), mean_logFC=round(t.logFC.mean(), 3),
            pct_positive=round(100 * (t.logFC > 0).mean(), 1), n_sig=int(sig.sum()),
            median_span_kb_sig=round(t.span_kb[sig].median(), 1),
            median_span_kb_all=round(t.span_kb.median(), 1)))
    pd.DataFrame(rows).to_csv(f"{RES}/snrna_gene_length_regression.tsv",
                              sep="\t", index=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--q1", default=f"{RES}/q1_raw_T2000_d25.tsv")
    ap.add_argument("--loo", default=f"{RES}/loo")
    args = ap.parse_args()

    meta = pd.read_csv(f"{DATA}/meta_T2000.tsv", sep="\t", index_col=0)
    full = pd.read_csv(args.q1, sep="\t").set_index("gene")
    hits = set(full.index[full["adj.P.Val"] < FDR])
    comps = list(meta[meta.cohort == "COMPARATOR"].index)

    LOO = {c: pd.read_csv(f"{args.loo}/q1_raw_drop_{c}.tsv", sep="\t").set_index("gene")
           for c in comps}

    # ---- 1. leave-one-section-out ------------------------------------------------
    rows = []
    for c in comps:
        t = LOO[c]
        j = full[["logFC"]].join(t[["logFC", "adj.P.Val"]], rsuffix="_loo").dropna()
        h = set(t.index[t["adj.P.Val"] < FDR])
        jh = j.loc[list(hits)]
        rows.append(dict(
            dropped=c, disease=meta.loc[c, "disease"], n_sig=len(h),
            n_up_CS=int((t.loc[list(h), "logFC"] > 0).sum()),
            r_all985=round(j.logFC.corr(j.logFC_loo), 4),
            r_192hits=round(jh.logFC.corr(jh.logFC_loo), 4),
            sign_concord_192=int((np.sign(jh.logFC) == np.sign(jh.logFC_loo)).sum()),
            retained_of_192=len(hits & h), new_hits=len(h - hits),
            median_abs_lfc_192=round(jh.logFC_loo.abs().median(), 3)))
    pd.DataFrame(rows).sort_values("n_sig").to_csv(
        f"{RES}/q1_leave_one_comparator_out.tsv", sep="\t", index=False)

    # ---- 1b. leave-one-disease-out -----------------------------------------------
    LODO, dropped_n = {}, {}
    for d in sorted(set(meta.loc[comps, "disease"])):
        secs = DISEASE_SECTIONS.get(d, [d])
        dropped_n[d] = len(secs)
        LODO[d] = (pd.read_csv(f"{args.loo}/q1_raw_dropdis_{d}.tsv", sep="\t").set_index("gene")
                   if d in DISEASE_SECTIONS else LOO[secs[0]])
    rows = []
    for d, t in LODO.items():
        j = full[["logFC"]].join(t[["logFC", "adj.P.Val"]], rsuffix="_loo").dropna()
        h = set(t.index[t["adj.P.Val"] < FDR])
        jh = j.loc[list(hits)]
        rows.append(dict(
            dropped_disease=d, n_sections_dropped=dropped_n[d],
            n_comparator_left=len(comps) - dropped_n[d], n_sig=len(h),
            n_up_CS=int((t.loc[list(h), "logFC"] > 0).sum()),
            r_all985=round(j.logFC.corr(j.logFC_loo), 4),
            r_192hits=round(jh.logFC.corr(jh.logFC_loo), 4),
            sign_concord_192=int((np.sign(jh.logFC) == np.sign(jh.logFC_loo)).sum()),
            retained_of_192=len(hits & h), new_hits=len(h - hits)))
    pd.DataFrame(rows).sort_values("n_sig").to_csv(
        f"{RES}/q1_leave_one_disease_out.tsv", sep="\t", index=False)

    # ---- 2. stable core + snRNA pooling ------------------------------------------
    core = set.intersection(*[set(LOO[c].index[LOO[c]["adj.P.Val"] < FDR]) for c in comps]) & hits
    core_dn = sorted(g for g in core if full.loc[g, "logFC"] < 0)
    core_up = sorted(g for g in core if full.loc[g, "logFC"] > 0)
    core_dis = set.intersection(*[set(t.index[t["adj.P.Val"] < FDR]) for t in LODO.values()]) & hits
    print(f"[core] section-stable {len(core)} ({len(core_dn)} down, {len(core_up)} up) "
          f"| disease-stable {len(core_dis)}")

    P = pool_snrna()

    G = pd.DataFrame(index=pd.Index(core_dn, name=None))
    G["logFC"] = full.loc[G.index, "logFC"]
    G["fdr"] = full.loc[G.index, "adj.P.Val"]
    G["worst_loo_fdr"] = pd.concat(
        [LOO[c].loc[G.index, "adj.P.Val"] for c in comps], axis=1).max(1)
    G["disease_stable"] = G.index.isin(core_dis)
    G["snrna_lfc_pool"] = P.reindex(G.index).lfc_pool
    G["snrna_fdr_pool"] = P.reindex(G.index).fdr_pool
    G["snrna_concordant"] = np.where(G.snrna_lfc_pool.isna(), "-",
                                     np.where(G.snrna_lfc_pool < 0, "yes", "no"))
    G.sort_values("logFC").to_csv(f"{RES}/q1_stable_down_115.tsv", sep="\t")

    # ---- 3. cross-platform concordance -------------------------------------------
    rng = np.random.default_rng(SEED)

    def decile_pools(genes):
        counts = P.loc[genes].expr_dec.value_counts()
        return counts, {d: P.index[P.expr_dec == d] for d in counts.index}

    def draw(counts, pools):
        return np.concatenate([rng.choice(pools[d], k, replace=False)
                               for d, k in counts.items()])

    def assess(genes, exp_sign, label):
        g = [x for x in genes if x in P.index]
        sub = P.loc[g]
        conc = int((np.sign(sub.lfc_pool) == exp_sign).sum())
        rep = int(((np.sign(sub.lfc_pool) == exp_sign) & (sub.fdr_pool < FDR)).sum())
        rho, prho = stats.spearmanr(full.loc[g, "logFC"], sub.lfc_pool)
        counts, pools = decile_pools(g)
        nc = np.array([(np.sign(P.loc[draw(counts, pools), "lfc_pool"]) == exp_sign).mean()
                       for _ in range(NPERM)])
        sign_ok = np.sign(P.lfc_pool) == exp_sign
        nr = np.array([(sign_ok.loc[(d := draw(counts, pools))]
                        & (P.loc[d, "fdr_pool"] < FDR)).mean() for _ in range(NPERM)])
        obs = conc / len(g)
        return dict(set=label, n_core=len(genes), n_in_snrna=len(g), concordant=conc,
                    frac_conc=round(obs, 3),
                    binom_p=stats.binomtest(conc, len(g), 0.5, alternative="greater").pvalue,
                    null_frac=round(nc.mean(), 3), null_p=round((nc >= obs).mean() + 1 / NPERM, 4),
                    mean_lfc_pool=round(sub.lfc_pool.mean(), 3),
                    ttest_p=stats.ttest_1samp(sub.lfc_pool, 0).pvalue,
                    rep_fdr5=rep, spearman_rho=round(rho, 3), spearman_p=prho,
                    conc_DCM=int((np.sign(sub.lfc_DCM) == exp_sign).sum()),
                    conc_ARVC=int((np.sign(sub.lfc_ARVC) == exp_sign).sum()),
                    null_rep_frac=round(nr.mean(), 3),
                    rep_frac=round(rep / len(g), 3),
                    rep_perm_p=round((nr >= rep / len(g)).mean() + 1 / NPERM, 4))

    A = pd.DataFrame([assess(core_dn, -1, f"stable down ({len(core_dn)})"),
                      assess(core_up, +1, f"stable up ({len(core_up)})")])
    A.to_csv(f"{RES}/q1_core_snrna_concordance.tsv", sep="\t", index=False)
    print(A[["set", "n_in_snrna", "concordant", "frac_conc", "null_frac",
             "null_p", "rep_fdr5", "spearman_rho"]].to_string(index=False))

    # ---- per-gene cross-platform table (figure panel A) ---------------------------
    cg = [g for g in core_dn + core_up if g in P.index]
    T = P.loc[cg, ["lfc_DCM", "fdr_DCM", "lfc_ARVC", "fdr_ARVC",
                   "lfc_pool", "se_pool", "fdr_pool", "AveExpr"]].copy()
    T.insert(0, "gene", T.index)
    T["spatial_logFC"] = full.loc[cg, "logFC"].values
    T["direction"] = np.where(T.spatial_logFC < 0, "down", "up")
    exp_sign = np.where(T.direction == "down", -1, 1)
    T["concordant"] = np.sign(T.lfc_pool) == exp_sign
    T["replicates_fdr5"] = T.concordant & (T.fdr_pool < FDR)
    T = T.sort_values(["direction", "spatial_logFC"], ascending=[True, True])
    T.to_csv(f"{RES}/q1_core_snrna_pergene.tsv", sep="\t", index=False)
    print(f"[cross-platform] {len(T)} of {len(core)} core genes testable in both arms")

    gene_length_scan()
    print("[cross-platform] wrote snrna_gene_length_regression.tsv "
          "(gene-length trend per snRNA-seq contrast)")


if __name__ == "__main__":
    main()
