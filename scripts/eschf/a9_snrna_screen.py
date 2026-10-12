#!/usr/bin/env python3
"""ESC HF figure v3, s.3 and A12: snRNA reliability screen. Run from repo root after a9_study_means.R.

Pre-specified (HANDOFF_figure_v3_ESCHF.md s.3), on the primary fit (all CM-dominant spots):
  Sets: (i) all FDR<5% spatial hits testable in snRNA (present in the pooled table);
        (ii) their QC-passing subset; (iii) the concordant set.
  Concordant set (s.5, as written): from ALL spatial hits, same sign of pooled snRNA logFC
        (CS vs DCM + ARVC, inverse-variance pooled) and pooled FDR < 5%. The v2 figure's set was
        built from the section-stable core only (q1_core_snrna_pergene.tsv, replicates_fdr5);
        it is carried alongside as in_concordant_core_v2.
  3a (applied): unreliable if |lfc_pool| > 3.
  3b (reported, not applied): source-specific near-absence if, for some snRNA study, the
        per-study mean log2 CPM (a9_study_mean_log2cpm.tsv.gz) is < 0 and the mean of the other
        four studies' per-study means is >= 2 higher.
  3c: descriptive sign agreement (lower / higher in CS; all hits and QC-passing) with the
        background fraction from the expression-decile-matched label-permutation null of
        a3a4_qc_nonfailing.py (NPERM = 2000, seed 0). Screened: 3a-failing genes removed from
        the hit set AND from the null pool. Unscreened values are recomputed with the same
        procedure and shown alongside the committed ones.
"""
import importlib.util
import numpy as np
import pandas as pd

V2 = "results/eschf/primary_v2"; OUT = "results/eschf"; FDR = 0.05; NPERM, SEED = 2000, 0; THR = 3.0
spec = importlib.util.spec_from_file_location("pool", "scripts/spatial_first/04_snrna_pool.py")
pool = importlib.util.module_from_spec(spec); spec.loader.exec_module(pool)
P = pool.pool_snrna()
P0 = pd.read_csv(f"{V2}/snrna_pooled_full.tsv", sep="\t").set_index("gene")
assert np.allclose(P.loc[P0.index, "lfc_pool"], P0.lfc_pool)
H = pd.read_csv(f"{V2}/q1_hits_annotated.tsv", sep="\t").set_index("gene")
H = H[H["adj.P.Val"] < FDR]
core = pd.read_csv(f"{V2}/q1_core_snrna_pergene.tsv", sep="\t").set_index("gene")
core_conc = set(core.index[core.replicates_fdr5])
SM = pd.read_csv(f"{OUT}/a9_study_mean_log2cpm.tsv.gz", sep="\t").set_index("gene")
STUD = list(SM.columns)
ch = {}
for tag in ("DCM", "HCM"):
    t = pd.read_csv(f"data/published/chaffin2022_cardiomyocyte_{tag}vsNF.tsv", sep="\t").set_index("gene")
    ch[tag] = t

# ---- per-gene table over hits testable in snRNA
T = H[["logFC", "adj.P.Val", "direction", "qc_pass", "in_qc_universe", "flag_locus_quant_outlier",
       "flag_segdup"]].rename(columns={"logFC": "spatial_logFC", "adj.P.Val": "spatial_fdr"})
T["testable_snrna"] = T.index.isin(P.index)
T = T.join(P[["lfc_pool", "fdr_pool", "lfc_DCM", "fdr_DCM", "lfc_ARVC", "fdr_ARVC"]])
for tag in ("DCM", "HCM"):
    T[f"chaffin_{tag}_logFC"] = ch[tag].logFC.reindex(T.index)
    T[f"chaffin_{tag}_bg_flag"] = ch[tag].bg_contamination_flag.reindex(T.index)
sgn = np.where(T.spatial_logFC < 0, -1, 1)
T["same_sign"] = np.sign(T.lfc_pool) == sgn
T["in_concordant_allhits"] = T.testable_snrna & T.same_sign & (T.fdr_pool < FDR)
T["in_concordant_core_v2"] = T.index.isin(core_conc)
T["fail_3a"] = T.testable_snrna & (T.lfc_pool.abs() > THR)
# 3b over the study means
def absence(g):
    if g not in SM.index:
        return pd.Series(dict(flag_3b=np.nan, absent_study=None, absent_mean=np.nan, others_mean=np.nan))
    m = SM.loc[g]; best = None
    for s in STUD:
        o = m.drop(s).mean()
        if m[s] < 0 and o - m[s] >= 2 and (best is None or m[s] < best[1]):
            best = (s, m[s], o)
    if best is None:
        s = m.idxmin(); return pd.Series(dict(flag_3b=False, absent_study=None, absent_mean=m[s], others_mean=m.drop(s).mean()))
    return pd.Series(dict(flag_3b=True, absent_study=best[0], absent_mean=best[1], others_mean=best[2]))
T = T.join(pd.DataFrame({g: absence(g) for g in T.index}).T)
T["flag_3b"] = T.flag_3b.astype("boolean")
T["in_concordant_allhits_screened"] = T.in_concordant_allhits & ~T.fail_3a
T["in_concordant_core_v2_screened"] = T.in_concordant_core_v2 & ~T.fail_3a
T = T.join(SM.add_prefix("mean_log2cpm_"))
T.sort_values(["direction", "spatial_logFC"]).to_csv(f"{OUT}/a9_snrna_screen.tsv", sep="\t")

tt = T[T.testable_snrna]
sets = {"all_hits_testable": tt, "qc_pass_hits_testable": tt[tt.qc_pass.astype(bool)],
        "concordant_allhits": tt[tt.in_concordant_allhits], "concordant_core_v2": tt[tt.in_concordant_core_v2]}
cnt = []
for k, d in sets.items():
    cnt.append(dict(set=k, n=len(d), n_fail_3a=int(d.fail_3a.sum()), genes_fail_3a=",".join(d.index[d.fail_3a]),
                    n_flag_3b=int(d.flag_3b.fillna(False).sum()), genes_flag_3b=",".join(d.index[d.flag_3b.fillna(False)]),
                    n_3a_and_3b=int((d.fail_3a & d.flag_3b.fillna(False)).sum()),
                    n_3a_only=int((d.fail_3a & ~d.flag_3b.fillna(False)).sum()),
                    n_3b_only=int((~d.fail_3a & d.flag_3b.fillna(False)).sum()),
                    n_3b_not_evaluable=int(d.flag_3b.isna().sum())))
C = pd.DataFrame(cnt)
# Chaffin genes beyond +/-3 (s.3a rationale check)
for tag in ("DCM", "HCM"):
    t = ch[tag]
    cnt.append(dict(set=f"chaffin_{tag}_vs_NF_published", n=len(t), n_fail_3a=int((t.logFC.abs() > THR).sum()),
                    genes_fail_3a=f"bg-flag-excluded: {int((t[t.bg_contamination_flag == 0].logFC.abs() > THR).sum())} of {int((t.bg_contamination_flag == 0).sum())}"))
pd.DataFrame(cnt).to_csv(f"{OUT}/a9_screen_counts.tsv", sep="\t", index=False)

# ---- 3b source-absence table (genes flagged among the three sets, plus FPGT sanity)
fl = T[T.testable_snrna & T.flag_3b.fillna(False)].copy()
fp = pd.DataFrame([absence("FPGT")], index=["FPGT"]).join(SM.add_prefix("mean_log2cpm_"))
assert bool(fp.flag_3b.iloc[0]), "FPGT not flagged by 3b"
A = pd.concat([fl[["direction", "spatial_logFC", "lfc_pool", "fail_3a", "in_concordant_allhits", "qc_pass",
                    "flag_3b", "absent_study", "absent_mean", "others_mean"] + [f"mean_log2cpm_{s}" for s in STUD]],
               fp.assign(note="sanity check (handoff s.3b reference gene; not a hit)")])
A.to_csv(f"{OUT}/a9_source_absence.tsv", sep="\t", index_label="gene")

# ---- suspect loci
loci = [("NDUFC2", "suspect"), ("KCTD14", "partner of NDUFC2"), ("CALM2", "suspect"), ("CALM1", "partner of CALM2"),
        ("CALM3", "partner of CALM2"), ("CALM1_CALM2_CALM3_sum", "summed CALM1+CALM2+CALM3 counts"),
        ("HSPB2", "suspect"), ("CRYAB", "partner of HSPB2")]
S = pd.DataFrame([dict(gene=g, role=r, in_snrna_counts=g in SM.index,
                       **({f"mean_log2cpm_{s}": SM.loc[g, s] for s in STUD} if g in SM.index else {}),
                       lfc_pool=P.lfc_pool.get(g, np.nan), spatial_logFC=H.logFC.get(g, np.nan),
                       spatial_hit=g in H.index) for g, r in loci])
S.to_csv(f"{OUT}/a9_suspect_loci.tsv", sep="\t", index=False)

# ---- 3c agreement, unscreened vs screened
def crossplat(genes, exp_sign, Pp, rng):
    g = [x for x in genes if x in Pp.index]; sub = Pp.loc[g]
    conc = int((np.sign(sub.lfc_pool) == exp_sign).sum())
    counts = sub.expr_dec.value_counts(); pools = {d: Pp.index[Pp.expr_dec == d] for d in counts.index}
    nc = np.array([(np.sign(Pp.loc[np.concatenate([rng.choice(pools[d], k, replace=False) for d, k in counts.items()]),
                                   "lfc_pool"]) == exp_sign).mean() for _ in range(NPERM)])
    return dict(n_testable=len(g), n_same_sign=conc, frac_same_sign=conc / len(g), background_frac=nc.mean(),
                n_concordant_fdr5=int(((np.sign(sub.lfc_pool) == exp_sign) & (sub.fdr_pool < FDR)).sum()))
Ps = P[P.lfc_pool.abs() <= THR]
rows = []
for screen, Pp in (("unscreened", P), ("screened_3a", Ps)):
    rng = np.random.default_rng(SEED)
    for direction, s in (("down", -1), ("up", 1)):
        for qc, d in (("all_hits", H), ("qc_passing_only", H[H.qc_pass.astype(bool)])):
            gs = list(d.index[d.direction == direction])
            rows.append(dict(screen=screen, direction=direction, gene_set=qc, **crossplat(gs, s, Pp, rng)))
R = pd.DataFrame(rows)
com = pd.read_csv(f"{V2}/a3a4_direction_descriptive.tsv", sep="\t")
R.to_csv(f"{OUT}/a9_concordance_screened.tsv", sep="\t", index=False)

# ---- A12 key genes
K = T.reindex(["JUP", "ACADVL", "IDH2"])[["spatial_logFC", "spatial_fdr", "lfc_pool", "fdr_pool", "lfc_DCM", "lfc_ARVC",
                                          "chaffin_DCM_logFC", "chaffin_HCM_logFC", "qc_pass", "in_qc_universe",
                                          "flag_locus_quant_outlier", "flag_segdup", "fail_3a", "flag_3b",
                                          "in_concordant_allhits", "in_concordant_allhits_screened",
                                          "in_concordant_core_v2", "in_concordant_core_v2_screened"]]
K.to_csv(f"{OUT}/a12_key_genes.tsv", sep="\t")
pd.set_option("display.width", 250)
print(C.to_string()); print(pd.DataFrame(cnt[-2:])[["set", "n", "n_fail_3a", "genes_fail_3a"]].to_string())
print(R.round(3).to_string()); print(com.head(8).to_string())
print(T.loc[T.fail_3a, ["direction", "spatial_logFC", "lfc_pool", "lfc_DCM", "lfc_ARVC", "chaffin_DCM_logFC", "chaffin_HCM_logFC", "qc_pass", "flag_3b", "absent_study"]].round(2).to_string())
print(K.T.to_string())
