#!/usr/bin/env python3
"""ESC HF revision: per-spot scores for the new primary (handoff_figure_fixes s.3-4).

Run from the repo root after a0_compare.py.  Inputs
  data/eschf/primary_v2/spotcounts_allcm_T2000.npz (+ _spots.tsv.gz, _genes.txt)
      thinned (2,000 UMI) counts for CS and comparator cardiomyocyte-dominant spots (a0_primary_export.py)
  data/eschf/eschf_perspot_T2000.tsv.gz   committed per-spot table (all CS + comparator spots)
Scores use the committed method: mean over member genes of log1p(CP10K) on the thinned counts.
Check first: the committed scores (old program list, ACADVL, IDH2, OXPHOS within the 985, CM-structural,
ECM7) are recomputed from the new count export and asserted equal to the committed table.
New scores: program = results/eschf/primary_v2/geneset_primary_program_qcpass_down.txt;
            OXPHOS  = KEGG hsa00190 within the new tested universe (same rule as before, new universe).
Output: data/eschf/primary_v2/eschf_perspot_T2000.tsv.gz -- the committed table with score_program and
score_oxphos_kegg (and their n_genes_) replaced on cardiomyocyte-dominant spots, NaN on other spots
(only cardiomyocyte-dominant spots enter A2), and data/eschf/primary_v2/a2_missing_genes.tsv.
"""
import numpy as np
import pandas as pd
import scipy.sparse as sp

D, V2, T, CP10K = "data/eschf/primary_v2", "results/eschf/primary_v2", 2000, 1e4


def read_set(path):
    return [l.strip() for l in open(path) if l.strip() and not l.startswith("#")]


X = sp.load_npz(f"{D}/spotcounts_allcm_T{T}.npz").tocsr()
spots = pd.read_csv(f"{D}/spotcounts_allcm_T{T}_spots.tsv.gz", sep="\t")
genes = read_set(f"{D}/spotcounts_allcm_T{T}_genes.txt")
assert X.shape == (len(spots), len(genes)), (X.shape, len(spots), len(genes))
assert (np.asarray(X.sum(1)).ravel() <= T).all(), "gene-subset counts exceed T"  # subset of a 2,000-UMI spot
gidx = {g: i for i, g in enumerate(genes)}
lc = X.multiply(CP10K / T).tocsr(); lc.data = np.log1p(lc.data)

old = pd.read_csv("data/eschf/eschf_perspot_T2000.tsv.gz", sep="\t")
key = ["section", "spot"]
assert not old.duplicated(key).any() and not spots.duplicated(key).any()
cm = old[old.cm_dominant].set_index(key)
assert set(map(tuple, spots[key].values)) == set(cm.index), "export spots != committed CM-dominant spots"
cm = cm.loc[list(map(tuple, spots[key].values))]


def score(gs, label, miss):
    present = [g for g in gs if g in gidx]
    miss += [dict(score=label, gene=g, reason="absent from exported gene space") for g in gs if g not in gidx]
    return np.asarray(lc[:, [gidx[g] for g in present]].mean(1)).ravel(), len(present)


u985 = set(read_set("metadata/tested_genes_985.txt"))
unew = set(read_set(f"{D}/universe_allcm.txt"))
kegg = read_set("metadata/geneset_kegg_oxphos_hsa00190.txt")
committed = {"program": read_set("metadata/geneset_primary_program_qcpass_down.txt"), "ACADVL": ["ACADVL"],
             "IDH2": ["IDH2"], "oxphos_kegg": [g for g in kegg if g in u985],
             "cm_structural": read_set("metadata/geneset_cmstructural.txt"), "ecm7": read_set("metadata/geneset_ecm7.txt")}
_m = []
for k, gs in committed.items():
    s, n = score(gs, k, _m)
    assert n == cm[f"n_genes_{k}"].iloc[0], (k, n)
    dmax = float(np.max(np.abs(s - cm[f"score_{k}"].values)))
    assert dmax < 1e-5, (k, dmax)
    print(f"[verify] committed score_{k}: recomputed, max |diff| = {dmax:.2e}")

miss = []
new = {"program": read_set(f"{V2}/geneset_primary_program_qcpass_down.txt"),
       "oxphos_kegg": [g for g in kegg if g in unew]}
out = old.copy()
for k, gs in new.items():
    s, n = score(gs, k, miss)
    m = pd.Series(s, index=pd.MultiIndex.from_frame(spots[key]))
    idx = pd.MultiIndex.from_frame(out[key])
    out[f"score_{k}"] = m.reindex(idx).values
    out[f"n_genes_{k}"] = n
    print(f"[new] score_{k}: {n} of {len(gs)} genes present")
assert out.loc[out.cm_dominant, "score_program"].notna().all()
pd.DataFrame(miss, columns=["score", "gene", "reason"]).to_csv(f"{D}/a2_missing_genes.tsv", sep="\t", index=False)
out.to_csv(f"{D}/eschf_perspot_T2000.tsv.gz", sep="\t", index=False)
print(len(new["oxphos_kegg"]), "OXPHOS genes in new universe;", len(committed["oxphos_kegg"]), "in the 985;",
      len(set(new["oxphos_kegg"]) & set(committed["oxphos_kegg"])), "shared;", len(miss), "missing")
