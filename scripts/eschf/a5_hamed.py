#!/usr/bin/env python3
"""ESC HF A5: comparison with Hamed et al. 2026 bulk-tissue proteomics. Run from repo root.

Different molecule (protein vs RNA), different resolution (bulk homogenate vs
cardiomyocyte-dominant spots), different cohort. Agreement and disagreement are reported
with equal prominence.

Inputs (data/external/): 44161_2026_880_MOESM4_ESM.xlsx (unimputed log2 matrix + metadata),
MOESM5 (per-etiology vs control DE), MOESM7 (NMF subgroup assignments).

Fixed rules (set before looking at any result):
  * sample names: matrix CT_n == metadata CI_n (cardiotoxicity); renamed CT_ -> CI_.
  * protein group -> gene: the Genes field is split on ';' and the row is mapped to each
    symbol. When a symbol maps to several rows, keep the row with the fewest member genes,
    then (MOESM4) the most quantified samples / (MOESM5) the smallest adjusted p.
  * per-etiology vs control (MOESM5): set genes present in the table; fraction with
    logFC < 0 vs background = tested spatial genes that are not FDR<5% hits; hypergeometric,
    both tails (lower tail = depletion, upper = enrichment).
  * Group 2 vs Group 1 (MOESM7 NMF_subgroup 2 vs 1), per protein OLS y ~ group + batch on
    the UNIMPUTED log2 values; fitted when >= 3 quantified samples per group; design rank
    asserted (rank-deficient proteins are set NA and counted); BH over fitted proteins.
    Age and sex are not available, so p-values are more liberal than the paper's.
  * cardiomyocyte content: per-sample mean z (z across the 122 grouped samples, missing
    ignored) of MYL2, MYL3, TNNT2, TNNI3, ACTC1, MYBPC3, ACTN2, TPM1; any protein group
    mapping to MYH9/MYH10/MYH11 is excluded. Refit y ~ group + batch + cm_content.
  * severity axis: Spearman of the CS-vs-comparator spatial logFC (985 tested genes) with
    the G2-G1 protein effect over overlapping genes; enrichment of CS-lower genes among
    G2-lower and G2-higher proteins (BH < 0.05), hypergeometric both tails against
    overlapping non-hit tested genes.
  * gene classes for colouring only: ECM = KEGG hsa04512 + the 7-gene ECM set + COL*;
    metabolic = KEGG hsa01100 + hsa00190 (ECM takes precedence).
The ECM half is NOT interpreted as biology: Hamed is bulk tissue including fibrotic regions,
whereas the spatial spots are cardiomyocyte-dominant by construction.
"""
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

EXT = "data/external"
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from eschf_paths import OUT, Q1, V2  # noqa: E402
FDR = 0.05
SARC = ["MYL2", "MYL3", "TNNT2", "TNNI3", "ACTC1", "MYBPC3", "ACTN2", "TPM1"]
PAN_MYO = {"MYH9", "MYH10", "MYH11"}
MIN_PER_GROUP = 3


def read_set(fn):
    return [l.strip() for l in open(f"metadata/{fn}") if l.strip() and not l.startswith("#")]


def explode(df, tiebreak_col, ascending):
    d = df.copy()
    d["gene_list"] = d.Genes.astype(str).str.split(";")
    d["n_members"] = d.gene_list.str.len()
    d = d.explode("gene_list").rename(columns={"gene_list": "gene"})
    d = d[d.gene.notna() & (d.gene != "nan") & (d.gene != "")]
    d = d.sort_values(["gene", "n_members", tiebreak_col], ascending=[True, True, ascending])
    return d.drop_duplicates("gene").set_index("gene")


def hyper(k, n, K, N):
    """k successes in n draws; K successes in population N. Returns (p_low, p_high)."""
    return stats.hypergeom.cdf(k, N, K, n), stats.hypergeom.sf(k - 1, N, K, n)


q1 = pd.read_csv(Q1, sep="\t").set_index("gene")
H = pd.read_csv(f"{OUT}/q1_hits_annotated.tsv", sep="\t").set_index("gene")
SETS = {"all": list(H.index[H.direction == "down"]),
        "qc_pass": list(H.index[(H.direction == "down") & H.qc_pass])}
bg_genes = [g for g in q1.index if g not in H.index]
rows_long, missing = [], []

# ------------------------------------------------------------------ 1. per-etiology vs control
x5 = pd.ExcelFile(f"{EXT}/44161_2026_880_MOESM5_ESM.xlsx")
vc = []
for et in x5.sheet_names:
    t = pd.read_excel(x5, et)
    t = t.rename(columns={f"logFC_{et}": "logFC", f"adj.P.Val_{et}": "adj"})
    t = explode(t, "adj", True)
    bg = [g for g in bg_genes if g in t.index]
    Kb = int((t.loc[bg, "logFC"] < 0).sum())
    for qc, gs in SETS.items():
        g = [x for x in gs if x in t.index]
        if et == x5.sheet_names[0]:
            missing += [dict(analysis="A5.1 Hamed MOESM5", gene_set=qc, gene=x) for x in gs if x not in t.index]
        k = int((t.loc[g, "logFC"] < 0).sum())
        plo, phi = hyper(k, len(g), k + Kb, len(g) + len(bg))
        r = dict(etiology=et, qc=qc, n_set=len(gs), n_in_proteome=len(g), n_lower_in_disease=k,
                 frac_lower=k / len(g), n_background=len(bg), bg_frac_lower=Kb / len(bg),
                 p_depletion=plo, p_enrichment=phi,
                 n_lower_fdr5=int(((t.loc[g, "logFC"] < 0) & (t.loc[g, "adj"] < FDR)).sum()),
                 bg_frac_lower_fdr5=float(((t.loc[bg, "logFC"] < 0) & (t.loc[bg, "adj"] < FDR)).mean()))
        vc.append(r)
        for st in ("n_in_proteome", "n_lower_in_disease", "frac_lower", "bg_frac_lower", "p_depletion", "p_enrichment"):
            rows_long.append(dict(analysis=f"hamed_{et}_vs_control", gene_universe="all_fdr5_hits",
                                  direction="down", qc=qc, statistic=st, value=r[st]))
pd.DataFrame(vc).to_csv(f"{OUT}/a5_hamed_vs_control.tsv", sep="\t", index=False)

# ------------------------------------------------------------------ 2. Group 2 vs Group 1
x4 = f"{EXT}/44161_2026_880_MOESM4_ESM.xlsx"
un = pd.read_excel(x4, "2) Unimputed")
meta = pd.read_excel(x4, "3) Metadata").set_index("Sample")
nmf = pd.read_excel(f"{EXT}/44161_2026_880_MOESM7_ESM.xlsx", "Metadata_NMF_Subgroups")
ren = lambda s: ("CI_" + s[3:]) if str(s).startswith("CT_") else s
un = un.rename(columns=ren)
nmf["Sample"] = nmf.Sample.map(ren)
assert set(nmf.Sample) <= set(un.columns) and set(nmf.Sample) <= set(meta.index)
assert nmf.NMF_subgroup.value_counts().to_dict() == {2: 63, 1: 59}
samples = list(nmf.Sample)
grp = (nmf.set_index("Sample").NMF_subgroup == 2).astype(float).loc[samples]
batch = (meta.loc[samples, "Batch number"] == "Batch II").astype(float)
un["n_quant"] = un[samples].notna().sum(1)
prot = explode(un, "n_quant", False)

# cardiomyocyte-content score
Y = prot[samples]
Z = Y.sub(Y.mean(1), axis=0).div(Y.std(1, ddof=1), axis=0)
sarc_rows = []
for g in SARC:
    cand = un[un.Genes.astype(str).str.split(";").apply(lambda L: g in L)]
    cand = cand[~cand.Genes.astype(str).str.split(";").apply(lambda L: bool(set(L) & PAN_MYO))]
    if not len(cand):
        missing.append(dict(analysis="A5.3 sarcomeric score", gene_set="sarcomeric", gene=g)); continue
    pick = cand.assign(nm=cand.Genes.astype(str).str.split(";").str.len()).sort_values(
        ["nm", "n_quant"], ascending=[True, False]).iloc[0]
    sarc_rows.append(dict(marker=g, protein_group=pick["Protein.Group"], genes=pick.Genes, n_quant=int(pick.n_quant)))
sarc_tab = pd.DataFrame(sarc_rows)
sarc_tab.to_csv(f"{OUT}/a5_sarcomeric_proteins_used.tsv", sep="\t", index=False)
SZ = un.set_index("Protein.Group").loc[sarc_tab.protein_group, samples]
SZ = SZ.sub(SZ.mean(1), axis=0).div(SZ.std(1, ddof=1), axis=0)
cm = SZ.mean(0, skipna=True)
assert cm.notna().all()
pd.DataFrame(dict(sample=samples, nmf_group=nmf.NMF_subgroup.values, batch=meta.loc[samples, "Batch number"].values,
                  etiology=meta.loc[samples, "Etiology"].values, cm_content=cm.values)).to_csv(
    f"{OUT}/a5_hamed_cm_content.tsv", sep="\t", index=False)


def fit_all(covars):
    res, n_def = [], 0
    for gene, yrow in Y.iterrows():
        y = yrow.values.astype(float)
        ok = ~np.isnan(y)
        n1, n2 = int((ok & (grp.values == 0)).sum()), int((ok & (grp.values == 1)).sum())
        if n1 < MIN_PER_GROUP or n2 < MIN_PER_GROUP:
            continue
        X = np.column_stack([np.ones(ok.sum()), grp.values[ok], batch.values[ok]] +
                            [c.values[ok] for c in covars])
        if np.linalg.matrix_rank(X) < X.shape[1]:
            n_def += 1; continue
        b, _, _, _ = np.linalg.lstsq(X, y[ok], rcond=None)
        resid = y[ok] - X @ b
        dfr = ok.sum() - X.shape[1]
        s2 = resid @ resid / dfr
        se = np.sqrt(s2 * np.linalg.inv(X.T @ X)[1, 1])
        tt = b[1] / se
        res.append(dict(gene=gene, effect=b[1], se=se, t=tt, p=2 * stats.t.sf(abs(tt), dfr), n_g1=n1, n_g2=n2))
    r = pd.DataFrame(res).set_index("gene")
    r["adj"] = multipletests(r.p, method="fdr_bh")[1]
    return r, n_def


raw, ndef_raw = fit_all([])
adj, ndef_adj = fit_all([cm])
print(f"[G2vsG1] fitted {len(raw)} proteins (rank-deficient {ndef_raw}); adjusted fitted {len(adj)} (rank-deficient {ndef_adj})")
sanity = {"POSTN": +1, "COL1A1": +1, "ASPN": +1, "CKM": -1, "TNNT2": -1}
san = pd.DataFrame([dict(gene=g, expected=("higher" if s > 0 else "lower"), effect=raw.loc[g, "effect"],
                         adj=raw.loc[g, "adj"], passes=bool(np.sign(raw.loc[g, "effect"]) == s))
                    for g, s in sanity.items()])
san.to_csv(f"{OUT}/a5_g2g1_sanity.tsv", sep="\t", index=False)
print(san.to_string(index=False))
assert san.passes.all(), "Group 2 vs Group 1 sanity check failed -- stop and report"

G = raw.add_suffix("_raw").join(adj.add_suffix("_cmadj"), how="outer")
G["n_members"] = prot.loc[G.index, "n_members"]
G["protein_group"] = prot.loc[G.index, "Protein.Group"]
ecm = set(read_set("geneset_kegg_ecm_receptor_hsa04512.txt")) | set(read_set("geneset_ecm7.txt"))
met = set(read_set("geneset_kegg_metabolic_hsa01100.txt")) | set(read_set("geneset_kegg_oxphos_hsa00190.txt"))
G["gene_class"] = ["ECM" if (g in ecm or g.startswith("COL")) else ("metabolic" if g in met else "other") for g in G.index]
G["st_logFC"] = q1.logFC.reindex(G.index)
G["st_hit"] = G.index.isin(H.index)
G["st_primary"] = G.index.isin(SETS["qc_pass"])
G.to_csv(f"{OUT}/a5_hamed_g2g1.tsv", sep="\t", index_label="gene")

# ------------------------------------------------------------------ 4. severity axis
sev = []
ov = [g for g in q1.index if g in raw.index]
ov_bg = [g for g in ov if g not in H.index]
for which, R in (("raw", raw), ("cm_adjusted", adj)):
    o = [g for g in ov if g in R.index]
    rho, p = stats.spearmanr(q1.loc[o, "logFC"], R.loc[o, "effect"])
    sev.append(dict(model=which, test="spearman_all_overlap", qc="all", n=len(o), value=rho, p=p))
    for qc, gs in SETS.items():
        g = [x for x in gs if x in R.index]
        bg = [x for x in ov_bg if x in R.index]
        if which == "raw":
            missing += [dict(analysis="A5.2/4 Hamed G2vsG1", gene_set=qc, gene=x) for x in gs if x not in R.index]
        rho_s, p_s = stats.spearmanr(q1.loc[g, "logFC"], R.loc[g, "effect"])
        sev.append(dict(model=which, test="spearman_cs_lower_set", qc=qc, n=len(g), value=rho_s, p=p_s))
        for half, cond in (("G2_lower", lambda d: (d.effect < 0) & (d.adj < FDR)),
                           ("G2_higher", lambda d: (d.effect > 0) & (d.adj < FDR))):
            k = int(cond(R.loc[g]).sum()); Kb = int(cond(R.loc[bg]).sum())
            plo, phi = hyper(k, len(g), k + Kb, len(g) + len(bg))
            sev.append(dict(model=which, test=f"{half}_count", qc=qc, n=len(g), value=k, p=np.nan,
                            frac=k / len(g), bg_frac=Kb / len(bg), n_bg=len(bg),
                            p_depletion=plo, p_enrichment=phi,
                            genes=",".join(sorted(np.array(g)[cond(R.loc[g]).values]))))
            if which == "raw":
                for st, v in ((f"{half}_n", k), (f"{half}_frac", k / len(g)), (f"{half}_bg_frac", Kb / len(bg)),
                              (f"{half}_p_enrichment", phi), (f"{half}_p_depletion", plo)):
                    rows_long.append(dict(analysis="hamed_G2_vs_G1", gene_universe="all_fdr5_hits", direction="down",
                                          qc=qc, statistic=st, value=v))
        if which == "raw":
            rows_long.append(dict(analysis="hamed_G2_vs_G1", gene_universe="all_fdr5_hits", direction="down",
                                  qc=qc, statistic="spearman_rho_set", value=rho_s))
# persistence of the G2-lower half after cardiomyocyte-content adjustment
for qc, gs in SETS.items():
    low = [g for g in gs if g in raw.index and raw.loc[g, "effect"] < 0 and raw.loc[g, "adj"] < FDR]
    keep = [g for g in low if g in adj.index and adj.loc[g, "effect"] < 0 and adj.loc[g, "adj"] < FDR]
    sev.append(dict(model="raw->cm_adjusted", test="G2_lower_persisting", qc=qc, n=len(low), value=len(keep),
                    genes=",".join(sorted(keep))))
    hi = [g for g in gs if g in raw.index and raw.loc[g, "effect"] > 0 and raw.loc[g, "adj"] < FDR]
    keep_hi = [g for g in hi if g in adj.index and adj.loc[g, "effect"] > 0 and adj.loc[g, "adj"] < FDR]
    sev.append(dict(model="raw->cm_adjusted", test="G2_higher_persisting", qc=qc, n=len(hi), value=len(keep_hi),
                    genes=",".join(sorted(keep_hi))))
S = pd.DataFrame(sev)
S.to_csv(f"{OUT}/a5_severity_axis.tsv", sep="\t", index=False)
print(S.drop(columns=["genes"], errors="ignore").to_string())
pd.DataFrame(missing).to_csv(f"{OUT}/a5_missing_genes.tsv", sep="\t", index=False)
rows_long.append(dict(analysis="hamed_G2_vs_G1", gene_universe="all_overlap", direction="all", qc="all",
                      statistic="spearman_rho_all_overlap", value=float(S.query("model=='raw' & test=='spearman_all_overlap'").value.iloc[0])))
rows_long.append(dict(analysis="hamed_G2_vs_G1", gene_universe="all_overlap", direction="all", qc="qc_pass",
                      statistic="spearman_rho_all_overlap",
                      value=stats.spearmanr(q1.loc[[g for g in ov if g in set(H.index[H.qc_pass]) or g not in H.index], "logFC"],
                                            raw.loc[[g for g in ov if g in set(H.index[H.qc_pass]) or g not in H.index], "effect"])[0]))

# ------------------------------------------------------------------ append to the A3 table
L = pd.read_csv(f"{OUT}/a3_qc_comparison_long.tsv", sep="\t")
L = pd.concat([L[~L.analysis.str.startswith("hamed")], pd.DataFrame(rows_long)], ignore_index=True)
L.to_csv(f"{OUT}/a3_qc_comparison_long.tsv", sep="\t", index=False)
W = L.pivot_table(index=["analysis", "gene_universe", "direction", "statistic"], columns="qc",
                  values="value", aggfunc="first", sort=False).reset_index()
W.columns.name = None
W.rename(columns={"all": "all_hits", "qc_pass": "qc_passing_only"}).to_csv(f"{OUT}/a3_qc_comparison.tsv", sep="\t", index=False)
print(pd.DataFrame(vc).to_string())
