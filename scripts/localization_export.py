#!/usr/bin/env python3
"""Per-spot program scores for the localization analysis (Foong CS Visium sections only).

Reproduces the thinning of spatial_first_export.py exactly: same inner-join gene space,
same per-section sc.pp.downsample_counts(random_state=0) at each T. Scores are per-spot
means of log1p(CP10K) over member genes -- no control-gene binning, no re-centering.
Region flags, dominant lineage and dist_gran_mm are NOT recomputed here; they are joined
locally from per_spot_regions.tsv.gz on (spot, section, T).
"""
import os, re, glob, json
import numpy as np, pandas as pd, scipy.sparse as sp
import scanpy as sc
from scipy.spatial import cKDTree

CFG = json.load(open("spatial_first_config.json"))
LOC = json.load(open("localization_config.json"))
RAW = "/gpfs/data/vangurilab/heartmap-cs/raw"
OUT = "out"; os.makedirs(OUT, exist_ok=True)
T_GRID = CFG["T_grid"]; PITCH = CFG["spot_pitch_um"]; RS = CFG["random_state"]
INNER = [l.strip() for l in open("inner_join_genes_genomewide.txt") if l.strip()]

SETS = {"program_qcpass": LOC["program_primary"]["genes"],
        "program_all": LOC["program_secondary"]["genes"],
        "cmstructural": LOC["cmstructural"]["genes"]}


def log(*a):
    print(*a, flush=True)


def collapse_symbols(X, symbols):
    sym = pd.Index(np.asarray(symbols, dtype=str)).str.upper()
    X = sp.csr_matrix(X)
    if not sym.has_duplicates:
        return X, sym
    codes, uniq = pd.factorize(sym)
    M = sp.csr_matrix((np.ones(len(codes)), (np.arange(len(codes)), codes)), shape=(len(codes), len(uniq)))
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
        X, sym = collapse_symbols(a.X, a.var_names)
        out.append(dict(section=sample, X=X, sym=sym, obs=np.asarray(a.obs_names, str), coords=coords))
        log(f"[load CS] {sample:12s} {X.shape[0]:5d} spots {X.shape[1]:6d} symbols coords={coords is not None}")
    return out


# ---- patient mapping: prefer a host-side sample_meta.tsv, else the section-number prefix
SM = f"{RAW}/sample_meta.tsv"
pat_map, pat_src = {}, "section_prefix"
if os.path.exists(SM):
    sm = pd.read_csv(SM, sep="\t", comment="#")
    if {"sample_id", "individual"} <= set(sm.columns):
        pat_map = dict(zip(sm.sample_id.astype(str), sm.individual.astype(str)))
        pat_src = "raw/sample_meta.tsv:individual"
log(f"[patient] source = {pat_src} | mapped ids = {len(pat_map)}")


def patient_of(sec):
    if sec in pat_map:
        return pat_map[sec]
    m = re.match(r"(CS[_-]?\d+)", sec)
    return m.group(1).replace("-", "_") if m else sec


S = load_cs()
for s in S:
    assert abs(s["X"].data - np.round(s["X"].data)).max() < 1e-6, f"{s['section']}: non-integer counts"
    miss = len(set(INNER) - set(s["sym"]))
    assert miss == 0, f"{s['section']}: {miss} inner-join genes absent -- gene space mismatch"

for s in S:
    s["X"] = sp.csr_matrix(s["X"][:, pd.Index(s["sym"]).get_indexer(INNER)])
    s["total"] = np.asarray(s["X"].sum(1)).ravel()

gidx = {g: i for i, g in enumerate(INNER)}
used = {k: [g for g in v if g in gidx] for k, v in SETS.items()}
json.dump({k: dict(n_requested=len(SETS[k]), n_used=len(used[k]), genes=used[k]) for k in SETS},
          open(f"{OUT}/localization_genes_used.json", "w"), indent=1)
log("[sets] " + " | ".join(f"{k}: {len(used[k])}/{len(SETS[k])}" for k in SETS))

tabs = []
for T in T_GRID:
    for s in S:
        m = s["total"] >= T
        n = int(m.sum())
        if n == 0:
            continue
        a = sc.AnnData(sp.csr_matrix(s["X"][m]).astype(np.float32))
        sc.pp.downsample_counts(a, counts_per_cell=T, random_state=RS)
        Xd = sp.csr_matrix(a.X)
        assert np.allclose(np.asarray(Xd.sum(1)).ravel(), T), f"{s['section']}: thinning missed T={T}"
        L = Xd.copy(); L.data = np.log1p(L.data / T * 1e4)          # log1p CP10K on thinned counts
        d = pd.DataFrame(dict(spot=s["obs"][m], section=s["section"], T=T,
                              patient=patient_of(s["section"]), umi_native=s["total"][m],
                              umi_thinned=float(T)))
        for k, gs in used.items():
            cols = [gidx[g] for g in gs]
            Z = np.asarray(L[:, cols].todense())
            d[f"score_{k}"] = Z.mean(1)
            d[f"ndet_{k}"] = (Z > 0).sum(1)
        if s["coords"] is not None:
            C = s["coords"][m]
            dd, _ = cKDTree(C).query(C, k=2)
            ppu = float(np.median(dd[:, 1])) / PITCH
            d["x_um"] = C[:, 0] / ppu; d["y_um"] = C[:, 1] / ppu; d["px_per_um"] = ppu
        else:
            d["x_um"] = np.nan; d["y_um"] = np.nan; d["px_per_um"] = np.nan
        tabs.append(d)
    log(f"[T={T}] sections done")

A = pd.concat(tabs, ignore_index=True)
A.to_csv(f"{OUT}/localization_perspot.tsv.gz", sep="\t", index=False)
summ = (A.groupby(["T", "section"]).agg(n_spots=("spot", "size"),
                                        patient=("patient", "first"),
                                        mean_program=("score_program_qcpass", "mean"),
                                        mean_cmstruct=("score_cmstructural", "mean")).reset_index())
summ.to_csv(f"{OUT}/localization_export_summary.tsv", sep="\t", index=False)
json.dump(dict(n_rows=int(len(A)), n_sections=int(A.section.nunique()), T_grid=T_GRID,
               patient_source=pat_src, n_patients=int(A.patient.nunique()),
               inner_join=len(INNER), random_state=RS),
          open(f"{OUT}/localization_export_meta.json", "w"), indent=1)
log("[done] rows=%d sections=%d patients=%d" % (len(A), A.section.nunique(), A.patient.nunique()))
