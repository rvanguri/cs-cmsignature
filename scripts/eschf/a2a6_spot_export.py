#!/usr/bin/env python3
"""ESC HF A2/A6 spot-level export (runs on the cluster next to the raw Visium data).

This is scripts/spatial_first/00_export_from_raw.py restricted to T = 2000, with the same
loaders, inner-join gene space, per-section thinning (sc.pp.downsample_counts,
random_state from spatial_first_config.json), pooled lineage z-scores over all 22
sections (CS + comparators) and the same CS-only immune-dominant
("granuloma") distance. Nothing about the region definition is changed.

Before writing anything it re-derives the committed discovery outputs and asserts they
match exactly:
  * per-section pseudobulk over SELECTED spots == data/spatial_first/pseudobulk_T2000.tsv.gz
  * n_cm_dominant / n_granulomatous / n_selected == data/spatial_first/region_selection.tsv
so the per-spot table below is provably the same spot universe as discovery.

Per-spot output (CS + comparator sections, all in-tissue spots retained at T=2000):
  dominant lineage, cm_dominant, selected, distance to nearest immune-dominant spot (mm,
  CS only), native UMI, pixel coordinates, and program scores. Program score = mean over
  member genes of log1p(CP10K) on the depth-matched (2,000 UMI) counts. No score_genes.

Usage (cluster):  python a2a6_spot_export.py <config.json> <comparators_dir> <genesets_dir> \
                     <pseudobulk_T2000.tsv.gz> <region_selection.tsv>
"""
import glob
import json
import os
import re
import shutil
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.spatial import cKDTree

CFG = json.load(open(sys.argv[1]))
COMP_DIR = sys.argv[2]
GS_DIR = sys.argv[3]
PB_REF = sys.argv[4]
RS_REF = sys.argv[5]
RAW = "/gpfs/data/vangurilab/heartmap-cs/raw"
OUT = "out"
os.makedirs(OUT, exist_ok=True)
os.makedirs(f"{OUT}/images", exist_ok=True)
T = 2000
assert T in CFG["T_grid"]
PITCH = CFG["spot_pitch_um"]
EXCL_MM = CFG["granuloma_exclusion_mm"]
LIN = CFG["lineages"]
GRAN_LIN = set(CFG["granuloma_lineages"])
RS = CFG["random_state"]
CP10K = 1e4


def log(*a):
    print(*a, flush=True)


def collapse_symbols(X, symbols):
    sym = pd.Index(np.asarray(symbols, dtype=str)).str.upper()
    X = sp.csr_matrix(X)
    if not sym.has_duplicates:
        return X, sym
    codes, uniq = pd.factorize(sym)
    M = sp.csr_matrix((np.ones(len(codes)), (np.arange(len(codes)), codes)),
                      shape=(len(codes), len(uniq)))
    return sp.csr_matrix(X @ M), pd.Index(uniq)


def _read_positions(files):
    for f in files:
        b = os.path.basename(f)
        if "tissue_positions" in b:
            p = pd.read_csv(f, compression="gzip" if f.endswith(".gz") else None)
            if "barcode" not in p.columns:
                p = pd.read_csv(f, header=None, compression="gzip" if f.endswith(".gz") else None,
                                names=["barcode", "in_tissue", "array_row", "array_col",
                                       "pxl_row_in_fullres", "pxl_col_in_fullres"])
            return p.set_index("barcode")
    return None


def load_cs():
    mats = sorted(glob.glob(f"{RAW}/foong_spatial/**/*filtered_feature_bc_matrix*.h5", recursive=True))
    allf = glob.glob(f"{RAW}/foong_spatial/**/*", recursive=True)
    out = []
    for m in mats:
        base = os.path.basename(m)
        gsm = (re.match(r"(GSM\d+)_", base) or [None, base.split("_")[0]])[1]
        sample = re.sub(r"\.h5$", "", re.sub(r"^GSM\d+_filtered_feature_bc_matrix[_-]?", "", base))
        sibs = [f for f in allf if os.path.dirname(f) == os.path.dirname(m)
                and os.path.basename(f).startswith(gsm + "_")]
        a = sc.read_10x_h5(m); a.var_names_make_unique()
        pos = _read_positions(sibs)
        coords = None
        if pos is not None:
            common = a.obs_names.intersection(pos.index)
            a = a[common].copy(); pos = pos.reindex(a.obs_names)
            if "in_tissue" in pos:
                keep = pos["in_tissue"].astype(int) == 1
                a = a[keep.values].copy(); pos = pos.reindex(a.obs_names)
            coords = pos[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy(dtype=float)
        # keep the H&E image + scale factors for the figure (panel A)
        for f in sibs:
            if re.search(r"(hires_image|scalefactors)", os.path.basename(f)):
                shutil.copy(f, f"{OUT}/images/{sample}__{os.path.basename(f)}")
        X, sym = collapse_symbols(a.X, a.var_names)
        out.append(dict(section=sample, cohort="CS", disease="CS", source="foong",
                        X=X, sym=sym, obs=np.asarray(a.obs_names, str), coords=coords))
        log(f"[load CS] {sample:12s} {X.shape[0]:5d} spots")
    return out


def load_comparators():
    out = []
    for d, dis in CFG["comparator_map"].items():
        sec = d.replace(" ", "_")
        p = f"{COMP_DIR}/{sec}"
        a = sc.read_10x_h5(f"{p}/filtered_feature_bc_matrix.h5"); a.var_names_make_unique()
        pos = _read_positions([f"{p}/spatial/tissue_positions.csv"])
        common = a.obs_names.intersection(pos.index)
        a = a[common].copy(); pos = pos.reindex(a.obs_names)
        keep = pos["in_tissue"].astype(int) == 1
        a = a[keep.values].copy(); pos = pos.reindex(a.obs_names)
        coords = pos[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy(dtype=float)
        X, sym = collapse_symbols(a.X, a.var_names)
        out.append(dict(section=sec, cohort="COMPARATOR", disease=dis, source="new_visium",
                        X=X, sym=sym, obs=np.asarray(a.obs_names, str), coords=coords))
        log(f"[load COMP] {sec:12s} {X.shape[0]:5d} spots")
    return out


def read_set(fn):
    return [l.strip().upper() for l in open(f"{GS_DIR}/{fn}") if l.strip() and not l.startswith("#")]


S = load_cs() + load_comparators()
inner = set.intersection(*[set(s["sym"]) for s in S])
inner = sorted(inner)
for s in S:
    s["X"] = sp.csr_matrix(s["X"][:, pd.Index(s["sym"]).get_indexer(inner)])
    s["total"] = np.asarray(s["X"].sum(1)).ravel()
gidx = {g: i for i, g in enumerate(inner)}
lin_names = list(LIN.keys())
lin_genes = {L: [g for g in gs if g in gidx] for L, gs in LIN.items()}

blocks, metas, keepmask = [], [], {}
for s in S:
    m = s["total"] >= T
    keepmask[s["section"]] = m
    a = ad.AnnData(sp.csr_matrix(s["X"][m]).astype(np.float32))
    sc.pp.downsample_counts(a, counts_per_cell=T, random_state=RS)
    Xd = sp.csr_matrix(a.X)
    assert np.allclose(np.asarray(Xd.sum(1)).ravel(), T)
    blocks.append(Xd)
    C = s["coords"][m] if s["coords"] is not None else np.full((int(m.sum()), 2), np.nan)
    metas.append(pd.DataFrame(dict(spot=s["obs"][m], section=s["section"], cohort=s["cohort"],
                                   disease=s["disease"], source=s["source"],
                                   umi_native=s["total"][m], px_x=C[:, 0], px_y=C[:, 1])))
Xall = sp.vstack(blocks).tocsr()
md = pd.concat(metas, ignore_index=True)

L1 = Xall.copy(); L1.data = np.log1p(L1.data)
scores = np.zeros((Xall.shape[0], len(lin_names)))
for j, L in enumerate(lin_names):
    Z = np.asarray(L1[:, [gidx[g] for g in lin_genes[L]]].todense())
    mu, sd = Z.mean(0), Z.std(0); sd[sd == 0] = 1.0
    scores[:, j] = ((Z - mu) / sd).mean(1)
dom = np.array(lin_names)[scores.argmax(1)]
md["dominant_lineage"] = dom
for j, L in enumerate(lin_names):
    md[f"z_{L}"] = scores[:, j]

md["dist_gran_mm"] = np.nan
md["px_per_um"] = np.nan
off = 0
for s in S:
    n = int(keepmask[s["section"]].sum()); sl = slice(off, off + n); off += n
    if s["coords"] is None or n < 10:
        continue
    C = s["coords"][keepmask[s["section"]]]
    d, _ = cKDTree(C).query(C, k=2)
    ppu = float(np.median(d[:, 1])) / PITCH
    md.loc[md.index[sl], "px_per_um"] = ppu
    if s["cohort"] != "CS":
        continue
    gm = np.isin(dom[sl], list(GRAN_LIN))
    if gm.sum():
        dg, _ = cKDTree(C[gm]).query(C, k=1)
        md.loc[md.index[sl], "dist_gran_mm"] = dg / ppu / 1000.0
md["cm_dominant"] = md.dominant_lineage.values == "Cardiomyocyte"
md["immune_dominant"] = md.dominant_lineage.isin(GRAN_LIN).values
md["fibroblast_dominant"] = md.dominant_lineage.values == "Fibroblast"
excl = (md.cohort.values == "CS") & (md.dist_gran_mm.values < EXCL_MM)
md["selected"] = md.cm_dominant.values & ~excl

# ---------------------------------------------------------------- verification vs discovery
ref_pb = pd.read_csv(PB_REF, sep="\t", index_col=0)
ref_rs = pd.read_csv(RS_REF, sep="\t"); ref_rs = ref_rs[ref_rs["T"] == T].set_index("section")
secs = [s["section"] for s in S]
pb = np.zeros((len(secs), len(inner)), dtype=np.int64)
for i, sec in enumerate(secs):
    m = (md.section.values == sec) & md.selected.values
    if m.sum():
        pb[i] = np.asarray(Xall[np.where(m)[0]].sum(0)).ravel().astype(np.int64)
pb = pd.DataFrame(pb, index=secs, columns=inner)
assert list(ref_pb.columns) == inner, "gene space differs from committed pseudobulk"
pb_eq = (pb.loc[ref_pb.index].values == ref_pb.values).all()
cnt = md.groupby("section").agg(n_cm_dominant=("cm_dominant", "sum"),
                                n_granulomatous=("immune_dominant", "sum"),
                                n_selected=("selected", "sum"))
cmp = cnt.join(ref_rs[["n_cm_dominant", "n_granulomatous", "n_selected"]], rsuffix="_ref")
cnt_eq = all((cmp[c] == cmp[f"{c}_ref"]).all() for c in ["n_cm_dominant", "n_granulomatous", "n_selected"])
cmp.to_csv(f"{OUT}/verify_region_counts.tsv", sep="\t")
log(f"[verify] pseudobulk identical: {pb_eq} | region counts identical: {cnt_eq}")
assert pb_eq and cnt_eq, "re-derived spot universe differs from committed discovery outputs"

# ---------------------------------------------------------------- program scores
sets = {
    "program": read_set("geneset_primary_program_qcpass_down.txt"),
    "ACADVL": ["ACADVL"], "IDH2": ["IDH2"],
    "oxphos_kegg": [g for g in read_set("geneset_kegg_oxphos_hsa00190.txt")
                    if g in set(read_set("tested_genes_985.txt"))],
    "cm_structural": read_set("geneset_cmstructural.txt"),
    "ecm7": read_set("geneset_ecm7.txt"),
}
keep_rows = np.isin(md.cohort.values, ["CS", "COMPARATOR"])
Xk = Xall[np.where(keep_rows)[0]]
lc = Xk.multiply(CP10K / T).tocsr(); lc.data = np.log1p(lc.data)
out = md[keep_rows].reset_index(drop=True)
miss = []
for k, gs in sets.items():
    present = [g for g in gs if g in gidx]
    for g in gs:
        if g not in gidx:
            miss.append(dict(score=k, gene=g, reason="absent from inner-join Visium gene space"))
    out[f"score_{k}"] = np.asarray(lc[:, [gidx[g] for g in present]].mean(1)).ravel()
    out[f"n_genes_{k}"] = len(present)
pd.DataFrame(miss, columns=["score", "gene", "reason"]).to_csv(f"{OUT}/a2_missing_genes.tsv", sep="\t", index=False)
cols = ["spot", "section", "cohort", "disease", "umi_native", "px_x", "px_y", "px_per_um",
        "dominant_lineage", "cm_dominant", "immune_dominant", "fibroblast_dominant", "selected",
        "dist_gran_mm"] + [f"z_{L}" for L in lin_names] + \
       [c for c in out.columns if c.startswith("score_") or c.startswith("n_genes_")]
out[cols].to_csv(f"{OUT}/eschf_perspot_T2000.tsv.gz", sep="\t", index=False)
json.dump(dict(T=T, n_spots=int(len(out)), inner_join=len(inner), verified_pseudobulk=bool(pb_eq),
               verified_region_counts=bool(cnt_eq), sets={k: len(v) for k, v in sets.items()}),
          open(f"{OUT}/eschf_export_meta.json", "w"), indent=1)
log("[done]")
