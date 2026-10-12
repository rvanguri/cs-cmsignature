#!/usr/bin/env python3
"""ESC HF A3 (QC-consistent statistics) and A4 (non-failing comparisons). Run from repo root.

Primary gene set = QC-passing hits (results/eschf/q1_hits_annotated.tsv, qc_pass == True).
Each statistic is computed twice, on all genes of a set and on its QC-passing subset:

  cross-platform (snRNA pooled CS vs DCM+ARVC; exactly the procedure of
  scripts/spatial_first/04_snrna_pool.py: inverse-variance pooling, expression-decile
  matched permutation null, NPERM = 2000, seed 0):
    * section-stable core (as published: stable down / stable up), and
    * all FDR<5% hits (down / up).
  Chaffin non-failing (scripts/spatial_first/06_published_reference.py, unchanged logic):
    CS-lower genes; reference rate = fraction of published CM differential genes LOWER in
    disease; one-sided binomial.
  Internal snRNA CS vs non-failing (data/snrna/snrna_CS_vs_NF.tsv, CS - NF from the same
    limma fit that produced CS vs DCM/ARVC): fraction of CS-lower genes lower in CS than NF;
    one-sided binomial against the fraction lower among all jointly testable genes
    (tested in the spatial contrast AND in the snRNA NF contrast).
Hamed rows are appended by a5_hamed.py.
"""
import importlib.util

import numpy as np
import pandas as pd
from scipy import stats

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from eschf_paths import OUT, RES, Q1, V2  # noqa: E402
FDR = 0.05
NPERM, SEED = 2000, 0

spec = importlib.util.spec_from_file_location("pool", "scripts/spatial_first/04_snrna_pool.py")
pool = importlib.util.module_from_spec(spec); spec.loader.exec_module(pool)

full = pd.read_csv(Q1, sep="\t").set_index("gene")
H = pd.read_csv(f"{OUT}/q1_hits_annotated.tsv", sep="\t").set_index("gene")
meta = pd.read_csv("data/spatial_first/meta_T2000.tsv", sep="\t", index_col=0)
comps = list(meta[meta.cohort == "COMPARATOR"].index)
hits = set(H.index)
qcp = set(H.index[H.qc_pass])
LOO = {c: pd.read_csv(f"{RES}/loo/q1_raw_drop_{c}.tsv", sep="\t").set_index("gene") for c in comps}
core = set.intersection(*[set(LOO[c].index[LOO[c]["adj.P.Val"] < FDR]) for c in comps]) & hits
core_dn = sorted(g for g in core if full.loc[g, "logFC"] < 0)
core_up = sorted(g for g in core if full.loc[g, "logFC"] > 0)
hit_dn = sorted(g for g in hits if full.loc[g, "logFC"] < 0)
hit_up = sorted(g for g in hits if full.loc[g, "logFC"] > 0)

P = pool.pool_snrna()
rng = np.random.default_rng(SEED)


def crossplat(genes, exp_sign):
    g = [x for x in genes if x in P.index]
    sub = P.loc[g]
    conc = int((np.sign(sub.lfc_pool) == exp_sign).sum())
    rep = int(((np.sign(sub.lfc_pool) == exp_sign) & (sub.fdr_pool < FDR)).sum())
    rho, prho = stats.spearmanr(full.loc[g, "logFC"], sub.lfc_pool)
    counts = P.loc[g].expr_dec.value_counts()
    pools = {d: P.index[P.expr_dec == d] for d in counts.index}
    nc = np.array([(np.sign(P.loc[np.concatenate([rng.choice(pools[d], k, replace=False)
                                                  for d, k in counts.items()]), "lfc_pool"])
                    == exp_sign).mean() for _ in range(NPERM)])
    obs = conc / len(g)
    return dict(n_set=len(genes), n_testable=len(g), n_concordant=conc, frac_concordant=obs,
                null_frac=nc.mean(), null_p=(nc >= obs).mean() + 1 / NPERM,
                binom_p_vs_half=stats.binomtest(conc, len(g), 0.5, alternative="greater").pvalue,
                n_replicate_fdr5=rep, spearman_rho=rho, spearman_p=prho)


rows = []
for uni, dn, up in (("section_stable_core", core_dn, core_up), ("all_fdr5_hits", hit_dn, hit_up)):
    for direction, genes, s in (("down", dn, -1), ("up", up, +1)):
        for qc_lab, gs in (("all", genes), ("qc_pass", [g for g in genes if g in qcp])):
            r = crossplat(gs, s)
            for k, v in r.items():
                rows.append(dict(analysis="cross_platform_snrna", gene_universe=uni,
                                 direction=direction, qc=qc_lab, statistic=k, value=v))

# ---- A4.1 Chaffin (same logic as 06_published_reference.py)
for tag in ("DCM", "HCM"):
    t = pd.read_csv(f"data/published/chaffin2022_cardiomyocyte_{tag}vsNF.tsv", sep="\t")
    t = t[t.bg_contamination_flag == 0].set_index("gene")
    base = float((t.logFC < 0).mean())
    for qc_lab, gs in (("all", hit_dn), ("qc_pass", [g for g in hit_dn if g in qcp])):
        g = [x for x in gs if x in t.index]
        k = int((t.loc[g, "logFC"] < 0).sum())
        for st, v in dict(n_set=len(gs), n_testable=len(g), n_lower_in_disease=k,
                          frac_lower_in_disease=k / len(g), ref_frac_lower_in_disease=base,
                          n_published_deg=len(t),
                          binom_p=stats.binomtest(k, len(g), base, alternative="greater").pvalue).items():
            rows.append(dict(analysis=f"chaffin_{tag}_vs_NF", gene_universe="all_fdr5_hits",
                             direction="down", qc=qc_lab, statistic=st, value=v))

# ---- A4.2 internal snRNA CS vs NF
nf = pd.read_csv("data/snrna/snrna_CS_vs_NF.tsv", sep="\t").set_index("gene")
joint = nf.index.intersection(full.index)
base_nf = float((nf.loc[joint, "logFC"] < 0).mean())
for qc_lab, gs in (("all", hit_dn), ("qc_pass", [g for g in hit_dn if g in qcp])):
    g = [x for x in gs if x in nf.index]
    k = int((nf.loc[g, "logFC"] < 0).sum())
    kf = int(((nf.loc[g, "logFC"] < 0) & (nf.loc[g, "adj.P.Val"] < FDR)).sum())
    for st, v in dict(n_set=len(gs), n_testable=len(g), n_lower_in_CS_vs_NF=k,
                      frac_lower_in_CS_vs_NF=k / len(g), n_jointly_testable=len(joint),
                      ref_frac_lower_jointly_testable=base_nf,
                      binom_p=stats.binomtest(k, len(g), base_nf, alternative="greater").pvalue,
                      n_lower_fdr5=kf).items():
        rows.append(dict(analysis="snrna_CS_vs_NF", gene_universe="all_fdr5_hits",
                         direction="down", qc=qc_lab, statistic=st, value=v))

long = pd.DataFrame(rows)
long.to_csv(f"{OUT}/a3_qc_comparison_long.tsv", sep="\t", index=False)
wide = long.pivot_table(index=["analysis", "gene_universe", "direction", "statistic"],
                        columns="qc", values="value", aggfunc="first", sort=False).reset_index()
wide.columns.name = None
wide = wide.rename(columns={"all": "all_hits", "qc_pass": "qc_passing_only"})
wide.to_csv(f"{OUT}/a3_qc_comparison.tsv", sep="\t", index=False)
wide[wide.analysis.str.startswith("chaffin")].to_csv(f"{OUT}/a4_published_nf.tsv", sep="\t", index=False)
wide[wide.analysis == "snrna_CS_vs_NF"].to_csv(f"{OUT}/a4_internal_cs_vs_nf.tsv", sep="\t", index=False)

# per-gene CS vs NF for the 18 concordant (replicating) genes
pg = pd.read_csv(f"{RES}/q1_core_snrna_pergene.tsv", sep="\t")
c18 = pg[pg.replicates_fdr5].set_index("gene")
print("[concordant set]", len(c18), "genes")
t18 = c18[["direction", "spatial_logFC", "lfc_pool", "fdr_pool"]].join(
    nf[["logFC", "adj.P.Val"]].rename(columns={"logFC": "lfc_CS_vs_NF", "adj.P.Val": "fdr_CS_vs_NF"}))
t18 = t18.join(H[["qc_pass", "in_qc_universe", "flag_locus_quant_outlier"]])
t18.to_csv(f"{OUT}/a4_snrna_nf_concordant18.tsv", sep="\t")

# ---- excluded genes by flag type (A3.4)
ex = H[~H.qc_pass].copy()


def why(r):
    if not r.in_qc_universe:
        return "not_evaluated_by_locus_qc"
    f = [n for n, c in (("locus_quant_outlier", "flag_locus_quant_outlier"), ("segdup", "flag_segdup"),
                        ("missing_annotation", "flag_missing_annotation"), ("readthrough", "flag_readthrough"),
                        ("source_specific_absence", "flag_source_specific_absence")) if r[c] is True or r[c] == "True"]
    return ";".join(f) if f else "other"


ex["exclusion_reason"] = ex.apply(why, axis=1)
ex["concordant_replicating_gene"] = ex.index.isin(c18.index)
ex = ex.reset_index()[["gene", "direction", "logFC", "adj.P.Val", "exclusion_reason", "in_qc_universe",
                       "flag_locus_quant_outlier", "flag_segdup", "flag_missing_annotation",
                       "flag_readthrough", "flag_source_specific_absence", "max_abs_length_dev_log2",
                       "concordant_replicating_gene"]].sort_values(["exclusion_reason", "logFC"])
ex.to_csv(f"{OUT}/a3_excluded_genes.tsv", sep="\t", index=False)
print(ex.exclusion_reason.value_counts().to_string())
print(wide.to_string())
print(t18.to_string())

# ---------------------------------------------------------------- one descriptive convention (handoff_figure_fixes s.6)
# Fraction of genes in the stated direction beside the background fraction for the four direction comparisons.
# p values are carried for all four (never selectively); none is used as a gate.
SPEC = [  # analysis, gene_universe, direction, n_testable, n_in_direction, frac, background frac, p
    ("cross_platform_snrna", "all_fdr5_hits", "down", "n_testable", "n_concordant", "frac_concordant", "null_frac", "null_p"),
    ("chaffin_DCM_vs_NF", "all_fdr5_hits", "down", "n_testable", "n_lower_in_disease", "frac_lower_in_disease", "ref_frac_lower_in_disease", "binom_p"),
    ("chaffin_HCM_vs_NF", "all_fdr5_hits", "down", "n_testable", "n_lower_in_disease", "frac_lower_in_disease", "ref_frac_lower_in_disease", "binom_p"),
    ("snrna_CS_vs_NF", "all_fdr5_hits", "down", "n_testable", "n_lower_in_CS_vs_NF", "frac_lower_in_CS_vs_NF", "ref_frac_lower_jointly_testable", "binom_p"),
]
STATED = {"cross_platform_snrna": "same sign as spatial (lower in CS by pooled snRNA)",
          "chaffin_DCM_vs_NF": "lower in DCM than NF", "chaffin_HCM_vs_NF": "lower in HCM than NF",
          "snrna_CS_vs_NF": "lower in CS than NF"}
drows = []
for an, uni, d, kn, kk, kf, kb, kp in SPEC:
    for col in ("all_hits", "qc_passing_only"):
        x = wide[(wide.analysis == an) & (wide.gene_universe == uni) & (wide.direction == d)].set_index("statistic")[col]
        drows.append(dict(comparison=an, gene_set=f"lower-in-CS hits ({col})", stated_direction=STATED[an],
                          n_testable=int(x[kn]), n_in_direction=int(x[kk]), frac_in_direction=x[kf],
                          background_frac=x[kb], p_value=x[kp],
                          p_method="label-permutation null" if an == "cross_platform_snrna" else "binomial vs background"))
pd.DataFrame(drows).to_csv(f"{OUT}/a3a4_direction_descriptive.tsv", sep="\t", index=False)
