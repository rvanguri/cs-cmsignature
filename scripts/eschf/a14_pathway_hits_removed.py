#!/usr/bin/env python3
"""ESC HF A14 (HANDOFF_pathway_hits_removed.md): a priori pathway scores with fit hits removed.
Run from the repo root after a10_pathway_levels.py. Specification fixed before running:

Gene sets. KEGG hsa00020 (TCA), hsa00071 (fatty acid degradation), hsa00190 (OXPHOS) from metadata/geneset_kegg_*.txt
  (the A10 files), intersected with the 982-gene tested universe. Versions:
    full            all in-universe members (as A10)
    minus_primary   members minus primary-fit hits (primary_v2/q1_hits_annotated.tsv; every row has FDR < 5%)
    minus_union     members minus primary hits and >=0.5 mm sensitivity-fit hits (results/eschf/q1_hits_annotated.tsv)
  A version with < 5 genes is reported and marked too_small_to_interpret.
Patient level (primary readout). Patient pseudobulk of all CM-dominant spots (8 CS, 6 comparator patients),
  voom-form log2 CPM (a10_patient_log2cpm_voomform.tsv). Score = mean over members. CS mean, comparator mean and
  range of patient scores; CS - comparator with Welch 95% CI and p. The full version is asserted equal to A10.
Expression-matched null (minus_primary version only). Expression = per-gene mean voom log2 CPM over all 14
  patients. Deciles = 10 quantile bins of that expression over the 982-gene universe. Pool = universe minus all
  primary hits minus all members of the pathway tested. 10,000 draws, numpy default_rng(0), pathways in the order
  TCA, fatty acid degradation, OXPHOS; each draw takes, for each decile in increasing order, as many pool genes of
  that decile (without replacement) as the minus_primary members have in it. Statistic = the same CS - comparator
  difference (equal to the mean over the set of per-gene CS - comparator differences).
  percentile = 100 * (#null < obs + 0.5 * #null == obs) / N;
  two-sided p = min(1, 2 * min((1 + #null <= obs) / (N + 1), (1 + #null >= obs) / (N + 1))).
Robustness (minus_primary version).
  Leave-one-gene-out: the patient-level difference dropping each member; range and the gene whose removal moves it most.
  Per-gene direction: CS - comparator mean log2 CPM per member; counts lower / higher in CS.
  By distance: per-spot score = mean log1p CP10K over members on the T = 2000 thinned CM-dominant spot counts (as
  A10); per populated bin, CS mean of patient means with 95% t-interval vs comparator mean and range of patient
  means. Per-patient OLS slope exactly as a2_perpatient.py: score ~ dist_mm + cm_structural + z(log10 UMI)
  [+ section indicators], patients with >= 2 populated bins, design rank asserted; mean slope with t 95% CI,
  sign counts, one-sample t p.
Holds (SUMMARY rule): Welch 95% interval excludes zero AND null percentile < 2.5 or > 97.5.
Outputs: results/eschf/a14_pathway_hits_removed.tsv, a14_pathway_null.tsv, a14_pathway_robustness.tsv,
  a14_pathway_null_draws.tsv.gz (every draw's genes, for the harness).
The figure is not changed.
"""
import numpy as np, pandas as pd, scipy.sparse as sp
from scipy import stats

D, V2, OUT, T, CP10K, NDRAW, MINSET = "data/eschf/primary_v2", "results/eschf/primary_v2", "results/eschf", 2000, 1e4, 10000, 5
def read_set(p): return [l.strip() for l in open(p) if l.strip() and not l.startswith("#")]
SETS = {"tca_cycle_hsa00020": "metadata/geneset_kegg_tca_hsa00020.txt",
        "fatty_acid_degradation_hsa00071": "metadata/geneset_kegg_fatty_acid_degradation_hsa00071.txt",
        "oxphos_hsa00190": "metadata/geneset_kegg_oxphos_hsa00190.txt"}
univ = read_set(f"{D}/universe_allcm.txt"); U = set(univ)
HP = pd.read_csv(f"{V2}/q1_hits_annotated.tsv", sep="\t"); assert (HP["adj.P.Val"] < 0.05).all()
HS = pd.read_csv(f"{OUT}/q1_hits_annotated.tsv", sep="\t"); assert (HS["adj.P.Val"] < 0.05).all()
P, S = set(HP.gene), set(HS.gene)

L = pd.read_csv(f"{OUT}/a10_patient_log2cpm_voomform.tsv", sep="\t").set_index("gene")
assert set(L.index) == U
mp = pd.read_csv(f"{V2}/section_patient_map.tsv", sep="\t"); arm = mp.drop_duplicates("patient").set_index("patient").arm
csp = [c for c in L.columns if arm[c] == "CS"]; cop = [c for c in L.columns if arm[c] == "COMPARATOR"]
assert len(csp) == 8 and len(cop) == 6
dg = L[csp].mean(1) - L[cop].mean(1)          # per-gene CS - comparator

def welch(gs):
    s = L.loc[list(gs)].mean(0); x, y = s[csp], s[cop]
    vx, vy = x.var(ddof=1) / len(x), y.var(ddof=1) / len(y); se = np.sqrt(vx + vy)
    df = (vx + vy) ** 2 / (vx ** 2 / (len(x) - 1) + vy ** 2 / (len(y) - 1)); d = x.mean() - y.mean(); q = stats.t.ppf(0.975, df)
    return dict(cs_mean=x.mean(), comparator_mean=y.mean(), comparator_min=y.min(), comparator_max=y.max(),
                diff_cs_minus_comparator=d, diff_ci_low=d - q * se, diff_ci_high=d + q * se, welch_df=df,
                welch_p=2 * stats.t.sf(abs(d / se), df))

A10 = pd.read_csv(f"{OUT}/a10_pathway_levels.tsv", sep="\t"); A10 = A10[A10.level == "patient_pseudobulk"].set_index("set")
rows, members = [], {}
for k, p in SETS.items():
    allm = read_set(p); full = [g for g in allm if g in U]
    vers = {"full": full, "minus_primary": [g for g in full if g not in P], "minus_union": [g for g in full if g not in P | S]}
    members[k] = dict(all=allm, **vers)
    for v, gs in vers.items():
        r = dict(set=k, version=v, n_genes=len(gs), too_small_to_interpret=len(gs) < MINSET,
                 removed_genes=",".join(sorted(set(full) - set(gs))), n_removed=len(set(full) - set(gs)),
                 members=",".join(gs), **welch(gs))
        rows.append(r)
    f = [r for r in rows if r["set"] == k and r["version"] == "full"][0]
    assert f["n_genes"] == A10.loc[k, "n_genes"] and np.isclose(f["diff_cs_minus_comparator"], A10.loc[k, "diff_cs_minus_comparator"]) \
        and np.isclose(f["welch_p"], A10.loc[k, "welch_p"]), f"full version != A10 for {k}"
R1 = pd.DataFrame(rows); R1.to_csv(f"{OUT}/a14_pathway_hits_removed.tsv", sep="\t", index=False)

# ---- expression-matched null
expr = L.mean(1); dec = pd.Series(pd.qcut(expr, 10, labels=False), index=L.index).astype(int)
rng = np.random.default_rng(0); nrows, drows = [], []
for k in SETS:
    gs = members[k]["minus_primary"]; obs = dg[gs].mean()
    pool = [g for g in univ if g not in P and g not in set(members[k]["all"])]
    need = dec[gs].value_counts().sort_index()
    byd = {d: np.array([g for g in pool if dec[g] == d]) for d in need.index}
    for d, n in need.items(): assert len(byd[d]) >= n, f"pool too small in decile {d} for {k}"
    null = np.empty(NDRAW)
    for i in range(NDRAW):
        pick = np.concatenate([rng.choice(byd[d], n, replace=False) for d, n in need.items()])
        null[i] = dg[pick].mean(); drows.append((k, i, ",".join(pick)))
    lo, hi, eq = (null < obs).sum(), (null > obs).sum(), (null == obs).sum()
    nrows.append(dict(set=k, version="minus_primary", n_genes=len(gs), n_pool=len(pool),
                      decile_counts=",".join(f"{d}:{n}" for d, n in need.items()), n_draws=NDRAW, seed=0,
                      observed_diff=obs, null_mean=null.mean(), null_sd=null.std(ddof=1),
                      null_q025=np.quantile(null, 0.025), null_q975=np.quantile(null, 0.975),
                      percentile=100 * (lo + 0.5 * eq) / NDRAW,
                      p_two_sided=min(1.0, 2 * min((1 + lo + eq) / (NDRAW + 1), (1 + hi + eq) / (NDRAW + 1)))))
R2 = pd.DataFrame(nrows); R2.to_csv(f"{OUT}/a14_pathway_null.tsv", sep="\t", index=False)
pd.DataFrame(drows, columns=["set", "draw", "genes"]).to_csv(f"{OUT}/a14_pathway_null_draws.tsv.gz", sep="\t", index=False)

# ---- robustness (minus_primary)
X = sp.load_npz(f"{D}/spotcounts_allcm_T{T}.npz").tocsr()
spots = pd.read_csv(f"{D}/spotcounts_allcm_T{T}_spots.tsv.gz", sep="\t")
genes = read_set(f"{D}/spotcounts_allcm_T{T}_genes.txt"); gidx = {g: i for i, g in enumerate(genes)}
lc = X.multiply(CP10K / T).tocsr(); lc.data = np.log1p(lc.data)
mi = pd.read_csv(f"{V2}/a2_model_input.tsv.gz", sep="\t")
mi["z_log10_umi"] = (mi.log10_umi - mi.log10_umi.mean()) / mi.log10_umi.std(ddof=1)
bins = pd.read_csv(f"{OUT}/a2_oxphos_level.tsv", sep="\t").bin.tolist()
rob = []
for k in SETS:
    gs = members[k]["minus_primary"]; base = welch(gs)["diff_cs_minus_comparator"]
    lo = {g: welch([h for h in gs if h != g])["diff_cs_minus_comparator"] for g in gs}
    gmax = max(lo, key=lambda g: abs(lo[g] - base))
    for g in gs:
        rob.append(dict(set=k, analysis="leave_one_gene_out", item=g, value=lo[g], change=lo[g] - base))
    rob.append(dict(set=k, analysis="leave_one_gene_out_summary", item="min|max|most_influential",
                    value=min(lo.values()), value2=max(lo.values()), gene=gmax, change=lo[gmax] - base))
    for g in gs:
        rob.append(dict(set=k, analysis="per_gene_cs_minus_comparator", item=g, value=dg[g]))
    rob.append(dict(set=k, analysis="per_gene_direction_counts", item="lower|higher in CS",
                    value=int((dg[gs] < 0).sum()), value2=int((dg[gs] > 0).sum())))
    sc = np.asarray(lc[:, [gidx[g] for g in gs]].mean(1)).ravel()
    Sx = spots[["spot", "section", "cohort"]].assign(score=sc)
    cs = mi.merge(Sx, on=["spot", "section"], how="left"); assert cs.score.notna().all() and len(cs) == len(mi)
    comp = Sx[Sx.cohort == "COMPARATOR"].groupby("section").score.mean(); assert len(comp) == 6
    pp = cs.groupby(["bin", "patient"]).score.mean().reset_index()
    for b in bins:
        v = pp[pp.bin == b].score; h = stats.t.ppf(0.975, len(v) - 1) * v.std(ddof=1) / np.sqrt(len(v))
        rob.append(dict(set=k, analysis="distance_bin", item=b, value=v.mean(), ci_low=v.mean() - h, ci_high=v.mean() + h,
                        n_patients=len(v), comparator_mean=comp.mean(), comparator_min=comp.min(), comparator_max=comp.max(),
                        cs_mean_below_comparator_min=v.mean() < comp.min(), cs_ci_high_below_comparator_min=v.mean() + h < comp.min()))
    nb = cs.groupby("patient").bin.nunique(); sl = []
    for p in sorted(nb.index[nb >= 2]):
        sub = cs[cs.patient == p]
        Xd = [np.ones(len(sub)), sub.dist_gran_mm.values, sub.score_cm_structural.values, sub.z_log10_umi.values]
        Xd += [(sub.section == s).values.astype(float) for s in sorted(sub.section.unique())[1:]]
        Xd = np.column_stack(Xd); assert np.linalg.matrix_rank(Xd) == Xd.shape[1], f"rank-deficient {p} {k}"
        b_ = np.linalg.lstsq(Xd, sub.score.values, rcond=None)[0][1]; sl.append(b_)
        rob.append(dict(set=k, analysis="patient_slope_per_mm", item=p, value=b_))
    sl = np.array(sl); h = stats.t.ppf(0.975, len(sl) - 1) * sl.std(ddof=1) / np.sqrt(len(sl))
    rob.append(dict(set=k, analysis="slope_summary", item="mean slope per mm", value=sl.mean(), ci_low=sl.mean() - h,
                    ci_high=sl.mean() + h, n_patients=len(sl), n_negative=int((sl < 0).sum()), n_positive=int((sl > 0).sum()),
                    t_p=stats.ttest_1samp(sl, 0).pvalue))
R3 = pd.DataFrame(rob); R3.to_csv(f"{OUT}/a14_pathway_robustness.tsv", sep="\t", index=False)

pd.set_option("display.width", 250)
print(R1.drop(columns=["members"]).round(4).to_string(index=False))
print(R2.round(4).to_string(index=False))
print(R3[R3.analysis.isin(["leave_one_gene_out_summary", "per_gene_direction_counts", "distance_bin", "slope_summary"])]
      .dropna(axis=1, how="all").round(4).to_string(index=False))
