#!/usr/bin/env python3
"""Old-vs-new comparison for the Kuppe-free rerun (primary_v3).

usage: v3_old_vs_new.py OLD_ROOT NEW_ROOT OUTDIR

OLD_ROOT is the committed tree, NEW_ROOT the v3 run tree (both with results/manuscript_results.tsv
written by scripts/manuscript_results.py). Writes
  OUTDIR/old_vs_new_manuscript_values.tsv   one row per harness claim (joined on claim text)
  OUTDIR/old_vs_new_gene_lists.tsv          one row per named gene list, with set differences
Flags are mechanical (no interpretation):
  changed                value differs (numeric: relative difference > 1e-9)
  direction_change       both numeric, non-zero, opposite sign
  significance_crossing  claim names a p-value/FDR and the two values fall on opposite sides of 0.05
  status_change          harness status differs (e.g. OK -> MISMATCH: manuscript value no longer reproduced)
"""
import glob
import os
import re
import sys

import numpy as np
import pandas as pd

OLD, NEW, OUT = sys.argv[1:4]
os.makedirs(OUT, exist_ok=True)


def num(x):
    try:
        v = float(str(x).replace("%", ""))
        return v if np.isfinite(v) else None
    except ValueError:
        return None


# ---- claims -------------------------------------------------------------------------------------
o = pd.read_csv(f"{OLD}/results/manuscript_results.tsv", sep="\t", dtype=str, keep_default_na=False)
n = pd.read_csv(f"{NEW}/results/manuscript_results.tsv", sep="\t", dtype=str, keep_default_na=False)
for d in (o, n):
    d["k"] = d.claim + "#" + d.groupby("claim").cumcount().astype(str)   # duplicate claim texts
m = o.merge(n, on="k", how="outer", suffixes=("_old", "_new"), indicator=True)
m["claim"] = m.claim_old.where(m.claim_old.notna(), m.claim_new)
PV = re.compile(r"(\bp\b|p-value|p value|\bp_|\bP\b|FDR|fdr|adj\.?P|q-value)")
rows = []
for _, r in m.iterrows():
    presence = {"both": "both", "left_only": "old only", "right_only": "v3 only"}[r._merge]
    vo, vn = r.get("recomputed_old"), r.get("recomputed_new")
    a, b = num(vo), num(vn)
    if presence != "both":
        changed = True
    elif a is not None and b is not None:
        changed = abs(a - b) > 1e-9 * max(1.0, abs(a), abs(b))
    else:
        changed = str(vo) != str(vn)
    dirch = bool(presence == "both" and a is not None and b is not None and a != 0 and b != 0
                 and np.sign(a) != np.sign(b))
    sigx = bool(presence == "both" and a is not None and b is not None and PV.search(r.claim)
                and 0 <= a <= 1 and 0 <= b <= 1 and ((a < 0.05) != (b < 0.05)))
    stch = presence == "both" and r.status_old != r.status_new
    rows.append(dict(claim=r.claim, presence=presence,
                     manuscript_text_value=r.get("manuscript_old") if presence != "v3 only" else "",
                     old_value=vo, v3_value=vn,
                     delta=(b - a) if (a is not None and b is not None) else "",
                     old_status=r.get("status_old"), v3_status=r.get("status_new"),
                     changed=changed, direction_change=dirch, significance_crossing=sigx,
                     status_change=stch, v3_note=r.get("note_new")))
C = pd.DataFrame(rows)
C.to_csv(f"{OUT}/old_vs_new_manuscript_values.tsv", sep="\t", index=False)


# ---- gene lists ---------------------------------------------------------------------------------
def rd(root, p):
    return pd.read_csv(f"{root}/{p}", sep="\t")


def tb(s):
    return s.astype(str).str.lower().isin(["true", "1"])


def stable_core(root, res):
    conc = rd(root, f"{res}/q1_core_snrna_concordance.tsv").set_index("set")
    k = int(re.search(r"\((\d+)\)", [i for i in conc.index if i.startswith("stable down")][0]).group(1))
    return set(rd(root, f"{res}/q1_stable_down_{k}.tsv").iloc[:, 0])


def readset(path):
    return {l.strip() for l in open(path) if l.strip() and not l.startswith("#")}


def lists(root):
    L = {}
    for arm, res, raw, hits in (("sensitivity", "results/spatial_first", "q1_raw_T2000_d25.tsv", "results/eschf/q1_hits_annotated.tsv"),
                                ("primary", "results/eschf/primary_v2", "q1_raw_allcm.tsv", "results/eschf/a0_hits_annotated.tsv")):
        q = rd(root, f"{res}/{raw}")
        L[f"{arm}: tested genes"] = set(q.gene)
        for dn, sgn in (("lower", -1), ("higher", 1)):
            L[f"{arm}: FDR<5% {dn} in CS"] = set(q.gene[(q["adj.P.Val"] < 0.05) & (np.sign(q.logFC) == sgn)])
        h = rd(root, hits)
        L[f"{arm}: FDR<5% hits passing locus QC"] = set(h.gene[tb(h.qc_pass)])
        L[f"{arm}: FDR<5% hits failing locus QC"] = set(h.gene[~tb(h.qc_pass)])
        L[f"{arm}: leave-one-out stable lower core"] = stable_core(root, res)
        pg = rd(root, f"{res}/q1_core_snrna_pergene.tsv")
        L[f"{arm}: core genes testable in snRNA"] = set(pg.gene)
        for d in ("down", "up"):
            L[f"{arm}: core concordant and snRNA FDR<5%, {d}"] = set(
                pg.gene[tb(pg.concordant) & tb(pg.replicates_fdr5) & (pg.direction == d)])
    L["sensitivity: QC-pass CS-lower program (metadata)"] = readset(f"{root}/metadata/geneset_primary_program_qcpass_down.txt")
    L["primary: QC-pass CS-lower program"] = readset(f"{root}/results/eschf/primary_v2/geneset_primary_program_qcpass_down.txt")
    s9 = rd(root, "results/eschf/a9_snrna_screen.tsv")
    for c in ("same_sign", "in_concordant_allhits", "in_concordant_core_v2", "fail_3a", "flag_3b",
              "in_concordant_allhits_screened", "in_concordant_core"):
        if c in s9.columns:
            L[f"a9 screen: {c}"] = set(s9.gene[tb(s9[c])])
    L["a3 excluded genes"] = set(rd(root, "results/eschf/a3_excluded_genes.tsv").gene)
    L["a12 key genes"] = set(rd(root, "results/eschf/a12_key_genes.tsv").gene)
    a14 = rd(root, "results/eschf/a14_pathway_hits_removed.tsv")
    for _, r in a14.iterrows():
        L[f"a14 removed genes: {r['set']} ({r['version']})"] = set(str(r.removed_genes).replace(";", ",").split(",")) - {"", "nan"}
    for f in ("data/eschf/a2_missing_genes.tsv", "results/eschf/a5_missing_genes.tsv", "results/eschf/a5b_missing_genes.tsv"):
        if os.path.exists(f"{root}/{f}"):
            d = rd(root, f)
            L[f"missing genes: {os.path.basename(f)}"] = set(d.iloc[:, 0].astype(str)) if len(d) else set()
    pv = rd(root, "results/eschf/a7_panel_values.tsv")
    q = pv.quantity[(pv.panel == "D") & (pv.quantity.astype(str).str.count(r"\|") == 2)].astype(str)  # GENE|series|stat
    L["a7 figure panel D genes"] = set(q.str.split("|").str[0])
    return L


Lo, Ln = lists(OLD), lists(NEW)
G = []
for k in list(dict.fromkeys(list(Lo) + list(Ln))):
    a, b = Lo.get(k), Ln.get(k)
    G.append(dict(gene_list=k, n_old="" if a is None else len(a), n_v3="" if b is None else len(b),
                  n_shared="" if a is None or b is None else len(a & b),
                  only_old="" if a is None or b is None else ",".join(sorted(a - b)),
                  only_v3="" if a is None or b is None else ",".join(sorted(b - a)),
                  identical=(a == b)))
pd.DataFrame(G).to_csv(f"{OUT}/old_vs_new_gene_lists.tsv", sep="\t", index=False)
print(len(C), "claims;", int(C.changed.sum()), "changed;", int(C.direction_change.sum()), "direction;",
      int(C.significance_crossing.sum()), "significance;", int(C.status_change.sum()), "status;",
      dict(C.presence.value_counts()))
print(len(G), "gene lists;", sum(not g["identical"] for g in G), "differ")
