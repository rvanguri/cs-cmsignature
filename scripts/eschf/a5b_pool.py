#!/usr/bin/env python3
"""ESC HF A5b (part 2): pool the transplant-only snRNA contrasts and test procurement.

Run after a5b_procurement.R. Pooling = 04_snrna_pool.py procedure (inverse-variance,
SE = |logFC / t|, BH; expression-decile matched null, NPERM = 2000, seed 0).
Sets: the 18 concordant genes (replicating at FDR<5% in the full-cohort pool), the
primary gene set (QC-passing CS-lower hits), and the QC-passing section-stable core halves.
VAD vs transplant: LVAD - TX difference of per-patient program score and ACADVL
(log1p CP10K), Welch 95% CI; myectomy and autopsy excluded. LVEF: Spearman with program
score over CS patients with a numeric LVEF (left-censored values excluded).
Disease duration is not in the Neyazi supplementary tables.
"""
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from eschf_paths import OUT, RES, Q1, V2  # noqa: E402
FDR, NPERM, SEED = 0.05, 2000, 0


def pool(path_fmt):
    sn = {}
    for tag in ("DCM", "ARVC"):
        t = pd.read_csv(path_fmt.format(tag), sep="\t").set_index("gene")
        t["se"] = (t.logFC / t.t).abs()
        sn[tag] = t
    P = pd.DataFrame({"lfc_DCM": sn["DCM"].logFC, "se_DCM": sn["DCM"].se, "lfc_ARVC": sn["ARVC"].logFC,
                      "se_ARVC": sn["ARVC"].se, "AveExpr": sn["DCM"].AveExpr}).dropna()
    wD, wA = 1 / P.se_DCM ** 2, 1 / P.se_ARVC ** 2
    P["lfc_pool"] = (wD * P.lfc_DCM + wA * P.lfc_ARVC) / (wD + wA)
    P["se_pool"] = np.sqrt(1 / (wD + wA))
    P["p_pool"] = 2 * stats.norm.sf((P.lfc_pool / P.se_pool).abs())
    P["fdr_pool"] = multipletests(P.p_pool, method="fdr_bh")[1]
    P["expr_dec"] = pd.qcut(P.AveExpr, 10, labels=False, duplicates="drop")
    return P


q1 = pd.read_csv(Q1, sep="\t").set_index("gene")
H = pd.read_csv(f"{OUT}/q1_hits_annotated.tsv", sep="\t").set_index("gene")
pg = pd.read_csv(f"{RES}/q1_core_snrna_pergene.tsv", sep="\t")
meta = pd.read_csv("data/spatial_first/meta_T2000.tsv", sep="\t", index_col=0)
comps = list(meta[meta.cohort == "COMPARATOR"].index)
LOO = {c: pd.read_csv(f"{RES}/loo/q1_raw_drop_{c}.tsv", sep="\t").set_index("gene") for c in comps}
core = set.intersection(*[set(LOO[c].index[LOO[c]["adj.P.Val"] < FDR]) for c in comps]) & set(H.index)
qcp = set(H.index[H.qc_pass])
SETS = {"concordant_core": sorted(pg.gene[pg.replicates_fdr5]),
        f"primary_qcpass_down_{int(((H.direction == 'down') & H.qc_pass).sum())}": sorted(H.index[(H.direction == "down") & H.qc_pass]),
        "stable_core_down_qcpass": sorted(g for g in core if g in qcp and q1.loc[g, "logFC"] < 0),
        "stable_core_up_qcpass": sorted(g for g in core if g in qcp and q1.loc[g, "logFC"] > 0)}

FULL = pool("data/snrna/snrna_CS_vs_{}.tsv")
TX = pool(OUT + "/a5b_snrna_txonly_CS_vs_{}.tsv")
FULL.to_csv(f"{OUT}/snrna_pooled_full.tsv", sep="\t", index_label="gene")
TX.to_csv(f"{OUT}/snrna_pooled_txonly.tsv", sep="\t", index_label="gene")
rows = []
for name, genes in SETS.items():
    for cohort, P in (("full", FULL), ("transplant_only", TX)):
        rng = np.random.default_rng(SEED)
        g = [x for x in genes if x in P.index]
        s = np.sign(q1.loc[g, "logFC"]).values
        conc = int((np.sign(P.loc[g, "lfc_pool"]).values == s).sum())
        rep = int(((np.sign(P.loc[g, "lfc_pool"]).values == s) & (P.loc[g, "fdr_pool"].values < FDR)).sum())
        exp_sign = int(np.sign(s.mean())) if abs(s.mean()) == 1 else None
        null_p = null_mean = np.nan
        if exp_sign is not None and name != "concordant_core":
            counts = P.loc[g].expr_dec.value_counts()
            pools = {d: P.index[P.expr_dec == d] for d in counts.index}
            nc = np.array([(np.sign(P.loc[np.concatenate([rng.choice(pools[d], k, replace=False)
                                                          for d, k in counts.items()]), "lfc_pool"]) == exp_sign).mean()
                           for _ in range(NPERM)])
            null_p = (nc >= conc / len(g)).mean() + 1 / NPERM
            null_mean = nc.mean()
        same = np.nan
        if cohort == "transplant_only":
            gg = [x for x in g if x in FULL.index]
            same = int((np.sign(FULL.loc[gg, "lfc_pool"]) == np.sign(TX.loc[gg, "lfc_pool"])).sum())
        rows.append(dict(analysis="A5b.1 transplant-only refit", gene_set=name, cohort=cohort, n_set=len(genes),
                         n_testable=len(g), n_concordant_with_spatial=conc, frac_concordant=conc / len(g),
                         null_p=null_p, null_mean=null_mean, n_replicate_fdr5=rep, n_same_sign_as_full=same,
                         spearman_full_vs_tx=(stats.spearmanr(FULL.loc[g, "lfc_pool"], TX.loc[g, "lfc_pool"])[0]
                                              if cohort == "transplant_only" else np.nan)))
pd.DataFrame(rows).to_csv(f"{OUT}/a5b_txonly_concordance.tsv", sep="\t", index=False)
t18 = FULL.loc[SETS["concordant_core"], ["lfc_pool", "fdr_pool"]].add_suffix("_full").join(
    TX.loc[SETS["concordant_core"], ["lfc_pool", "fdr_pool"]].add_suffix("_txonly"))
t18.insert(0, "spatial_logFC", q1.loc[t18.index, "logFC"])
t18.to_csv(f"{OUT}/a5b_concordant18_txonly.tsv", sep="\t", index_label="gene")

# VAD vs transplant within CS; LVEF association
S = pd.read_csv(f"{OUT}/a5b_cs_patient_scores.tsv", sep="\t")
proc = []
for v in ("program_score", "ACADVL"):
    a, b = S.loc[S.procurement == "LVAD", v], S.loc[S.procurement == "TX", v]
    va, vb = a.var(ddof=1) / len(a), b.var(ddof=1) / len(b)
    se = np.sqrt(va + vb)
    dfw = (va + vb) ** 2 / (va ** 2 / (len(a) - 1) + vb ** 2 / (len(b) - 1))
    d = a.mean() - b.mean()
    q = stats.t.ppf(0.975, dfw)
    proc.append(dict(analysis="A5b.2 LVAD vs TX within CS", variable=v, n_LVAD=len(a), n_TX=len(b),
                     mean_LVAD=a.mean(), mean_TX=b.mean(), diff_LVAD_minus_TX=d, ci_low=d - q * se,
                     ci_high=d + q * se, welch_df=dfw, p=2 * stats.t.sf(abs(d / se), dfw),
                     mannwhitney_p=stats.mannwhitneyu(a, b, alternative="two-sided").pvalue))
L = S.dropna(subset=["lvef"])
for v in ("program_score", "ACADVL"):
    r, p = stats.spearmanr(L.lvef, L[v])
    proc.append(dict(analysis="A5b.3 LVEF association (all CS, numeric LVEF)", variable=v, n=len(L),
                     n_excluded_censored_or_missing=int(S.lvef.isna().sum()), spearman_rho=r, p=p))
proc.append(dict(analysis="A5b.3 disease duration", variable="disease_duration",
                 note="not available in Neyazi supplementary tables S1/S2"))
pd.concat([pd.DataFrame(rows), pd.DataFrame(proc)], ignore_index=True).to_csv(
    f"{OUT}/a5b_procurement.tsv", sep="\t", index=False)
print(pd.DataFrame(rows).to_string()); print(pd.DataFrame(proc).to_string()); print(t18.round(3).to_string())
