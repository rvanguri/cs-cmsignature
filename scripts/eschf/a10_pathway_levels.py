#!/usr/bin/env python3
"""ESC HF figure v3, A10: a priori pathway levels (KEGG fatty acid degradation hsa00071, TCA cycle hsa00020;
OXPHOS hsa00190 through the identical calculation for comparison). Run from repo root after a10_patient_log2cpm.R.

Pre-specified:
  Gene sets = KEGG members (metadata/geneset_kegg_*.txt, rest.kegg.jp) intersected with the 982-gene tested
  universe of the primary fit.
  Per-spot score = mean over set genes of log1p(CP10K) on the thinned (T = 2000) CM-dominant spot counts
  (identical to score_oxphos_kegg, a0_perspot_scores.py; OXPHOS recomputed and asserted equal).
  Distance bins and comparator reference exactly as a2a6_analysis.py / a2_perpatient.py (OXPHOS check):
  per populated bin, CS mean of patient means with 95% t-interval; comparator mean and range of patient means
  (all CM-dominant spots); flag whether the CS bin mean / upper CI falls below the comparator minimum.
  Patient level: per-patient pseudobulk log2 CPM (voom form, TMM; a10_patient_log2cpm_voomform.tsv), averaged
  over set genes; CS - comparator difference with Welch 95% CI (8 vs 6 patients).
Panel B is not changed by this analysis.
"""
import numpy as np, pandas as pd, scipy.sparse as sp
from scipy import stats

D, V2, OUT, T, CP10K = "data/eschf/primary_v2", "results/eschf/primary_v2", "results/eschf", 2000, 1e4
def read_set(p): return [l.strip() for l in open(p) if l.strip() and not l.startswith("#")]
X = sp.load_npz(f"{D}/spotcounts_allcm_T{T}.npz").tocsr()
spots = pd.read_csv(f"{D}/spotcounts_allcm_T{T}_spots.tsv.gz", sep="\t")
genes = read_set(f"{D}/spotcounts_allcm_T{T}_genes.txt"); gidx = {g: i for i, g in enumerate(genes)}
lc = X.multiply(CP10K / T).tocsr(); lc.data = np.log1p(lc.data)
univ = set(read_set(f"{D}/universe_allcm.txt"))
H = pd.read_csv(f"{V2}/q1_hits_annotated.tsv", sep="\t").set_index("gene"); H = H[H["adj.P.Val"] < 0.05]
SETS = {"fatty_acid_degradation_hsa00071": "metadata/geneset_kegg_fatty_acid_degradation_hsa00071.txt",
        "tca_cycle_hsa00020": "metadata/geneset_kegg_tca_hsa00020.txt",
        "oxphos_hsa00190": "metadata/geneset_kegg_oxphos_hsa00190.txt"}
mem, sc = [], {}
for k, p in SETS.items():
    gs = read_set(p); inu = [g for g in gs if g in univ]
    assert all(g in gidx for g in inu)
    sc[k] = np.asarray(lc[:, [gidx[g] for g in inu]].mean(1)).ravel()
    for g in gs:
        mem.append(dict(set=k, gene=g, in_universe=g in univ, primary_hit=g in H.index,
                        hit_direction=H.direction.get(g), hit_logFC=H.logFC.get(g)))
M = pd.DataFrame(mem); M.to_csv(f"{OUT}/a10_pathway_members.tsv", sep="\t", index=False)
S = spots[["spot", "section", "cohort"]].copy()
for k in SETS: S[k] = sc[k]
mi = pd.read_csv(f"{V2}/a2_model_input.tsv.gz", sep="\t")
cs = mi.merge(S, on=["spot", "section"], how="left"); assert cs[list(SETS)].notna().all().all()
assert np.abs(cs.oxphos_hsa00190 - cs.score_oxphos_kegg).max() < 1e-5   # identical to the committed OXPHOS score
cr = pd.read_csv(f"{V2}/a2_comparator_reference_perpatient.tsv", sep="\t").set_index("patient")
comp = S[S.cohort == "COMPARATOR"].groupby("section")[list(SETS)].mean()   # comparator patient = section
assert set(comp.index) == set(cr.index) and np.allclose(comp.loc[cr.index, "oxphos_hsa00190"], cr.score_oxphos_kegg, atol=1e-6)
ox0 = pd.read_csv(f"{OUT}/a2_oxphos_level.tsv", sep="\t")
rows = []
for k in SETS:
    pp = cs.groupby(["bin", "patient"])[k].mean().reset_index()
    for b in ox0.bin:
        v = pp[pp.bin == b][k]; h = stats.t.ppf(0.975, len(v) - 1) * v.std(ddof=1) / np.sqrt(len(v))
        c = comp[k]
        rows.append(dict(level="distance_bin", set=k, bin=b, n_patients=len(v), n_spots=int((cs.bin == b).sum()),
                         n_genes=int(M[(M.set == k)].in_universe.sum()), cs_mean_of_patient_means=v.mean(),
                         cs_ci_low=v.mean() - h, cs_ci_high=v.mean() + h, comparator_mean=c.mean(),
                         comparator_min=c.min(), comparator_max=c.max(), n_comparator_patients=len(c),
                         cs_mean_below_comparator_min=v.mean() < c.min(), cs_ci_high_below_comparator_min=v.mean() + h < c.min()))
r = pd.DataFrame(rows)
o = r[r.set == "oxphos_hsa00190"].set_index("bin")
assert np.allclose(o.loc[ox0.bin, "cs_mean_of_patient_means"], ox0.cs_mean_of_patient_means) and \
       np.allclose(o.loc[ox0.bin, "comparator_min"], ox0.comparator_min)
# patient-level pseudobulk
L = pd.read_csv(f"{OUT}/a10_patient_log2cpm_voomform.tsv", sep="\t").set_index("gene")
mp = pd.read_csv(f"{V2}/section_patient_map.tsv", sep="\t"); arm = mp.drop_duplicates("patient").set_index("patient").arm
for k in SETS:
    gs = M[(M.set == k) & M.in_universe].gene; s = L.loc[gs].mean(0)
    x, y = s[arm.reindex(s.index) == "CS"], s[arm.reindex(s.index) == "COMPARATOR"]
    vx, vy = x.var(ddof=1) / len(x), y.var(ddof=1) / len(y); se = np.sqrt(vx + vy)
    df = (vx + vy) ** 2 / (vx ** 2 / (len(x) - 1) + vy ** 2 / (len(y) - 1)); d = x.mean() - y.mean(); q = stats.t.ppf(0.975, df)
    rows.append(dict(level="patient_pseudobulk", set=k, bin="all CM-dominant spots", n_patients=len(x),
                     n_comparator_patients=len(y), n_genes=len(gs), cs_mean_of_patient_means=x.mean(),
                     comparator_mean=y.mean(), comparator_min=y.min(), comparator_max=y.max(),
                     diff_cs_minus_comparator=d, diff_ci_low=d - q * se, diff_ci_high=d + q * se,
                     welch_p=2 * stats.t.sf(abs(d / se), df)))
R = pd.DataFrame(rows); R.to_csv(f"{OUT}/a10_pathway_levels.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
print(M.groupby("set").agg(n=("gene", "size"), in_u=("in_universe", "sum"), hits=("primary_hit", "sum")).to_string())
print(M[M.primary_hit][["set", "gene", "hit_direction", "hit_logFC"]].round(2).to_string(index=False))
print(R.drop(columns=["n_spots"]).round(3).to_string(index=False))
