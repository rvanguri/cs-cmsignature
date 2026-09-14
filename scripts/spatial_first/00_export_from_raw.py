#!/usr/bin/env python3
"""Spatial-first export: genome-wide, depth-matched, section-level pseudobulk over
cardiomyocyte-dominant spots for 21 Visium sections (12 Foong CS, 6 comparators, 4 Kuppe normal).

Every threshold is read from spatial_first_config.json, fixed before any gene names are seen.
Region definition is COMPOSITIONAL (cardiomyocyte-argmax), not histologic; comparator sections
carry NO lesion exclusion.
"""
import os, re, sys, glob, json, gzip
import numpy as np, pandas as pd, scipy.sparse as sp
import scanpy as sc, anndata as ad
from scipy.spatial import cKDTree

CFG = json.load(open("spatial_first_config.json"))
RAW = ""
COMP_DIR = "comparators"
OUT = "out"; os.makedirs(OUT, exist_ok=True)
T_GRID = CFG["T_grid"]; PITCH = CFG["spot_pitch_um"]; EXCL_MM = CFG["granuloma_exclusion_mm"]
LIN = CFG["lineages"]; GRAN_LIN = set(CFG["granuloma_lineages"])
MIN_SPOTS_EVALUABLE = 50   # pre-specified: a section with fewer retained spots is not evaluable
RS = CFG["random_state"]


def log(*a):
    print(*a, flush=True)


def collapse_symbols(X, symbols):
    """Sum counts across duplicated gene symbols -> unique uppercase symbol space."""
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
            if "barcode" not in p.columns:   # headerless Space Ranger v1
                p = pd.read_csv(f, header=None, compression="gzip" if f.endswith(".gz") else None,
                                names=["barcode", "in_tissue", "array_row", "array_col",
                                       "pxl_row_in_fullres", "pxl_col_in_fullres"])
            return p.set_index("barcode")
    return None


def load_cs():
    """GEO flat Space Ranger files: group siblings by shared GSM prefix."""
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
        out.append(dict(section=sample, cohort="CS", disease="CS", source="foong",
                        X=X, sym=sym, obs=np.asarray(a.obs_names, str), coords=coords))
        log(f"[load CS] {sample:12s} {X.shape[0]:5d} spots {X.shape[1]:6d} symbols coords={coords is not None}")
    return out


def load_normal():
    out = []
    for f in sorted(glob.glob(f"{RAW}/kuppe_control/*.h5ad")):
        a = ad.read_h5ad(f)
        sec = "NORMAL_" + os.path.basename(f)[:8]
        vn = a.var_names.astype(str)
        e2s = None
        if pd.Index(vn).str.startswith("ENSG").mean() > 0.3:
            for col in ("feature_name", "gene_symbols", "gene_symbol", "gene_name", "Symbol", "symbol"):
                if col in a.var.columns:
                    e2s = dict(zip(vn, a.var[col].astype(str))); break
        src = a.raw if a.raw is not None else a
        rsym = np.asarray(src.var_names, str)
        if "feature_name" in getattr(src, "var", pd.DataFrame()).columns:
            rsym = src.var["feature_name"].astype(str).values
        elif e2s is not None:
            rsym = np.array([e2s.get(s, s) for s in rsym])
        keep = (a.obs["in_tissue"].astype(str) == "1").values if "in_tissue" in a.obs.columns \
            else np.ones(a.n_obs, bool)
        X = sp.csr_matrix(src.X)[keep]
        X, sym = collapse_symbols(X, rsym)
        coords = np.asarray(a.obsm["spatial"], float)[keep] if "spatial" in a.obsm else None
        out.append(dict(section=sec, cohort="NORMAL", disease="normal_donor", source="kuppe",
                        X=X, sym=sym, obs=np.asarray(a.obs_names, str)[keep], coords=coords))
        log(f"[load NORMAL] {sec:20s} {X.shape[0]:5d} spots {X.shape[1]:6d} symbols coords={coords is not None}")
    return out


def load_comparators():
    out = []
    for d, dis in CFG["comparator_map"].items():
        sec = d.replace(" ", "_")
        p = f"{COMP_DIR}/{sec}"
        a = sc.read_10x_h5(f"{p}/filtered_feature_bc_matrix.h5"); a.var_names_make_unique()
        pos = _read_positions([f"{p}/spatial/tissue_positions.csv"])
        coords = None
        if pos is not None:
            common = a.obs_names.intersection(pos.index)
            a = a[common].copy(); pos = pos.reindex(a.obs_names)
            keep = pos["in_tissue"].astype(int) == 1
            a = a[keep.values].copy(); pos = pos.reindex(a.obs_names)
            coords = pos[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy(dtype=float)
        X, sym = collapse_symbols(a.X, a.var_names)
        out.append(dict(section=sec, cohort="COMPARATOR", disease=dis, source="new_visium",
                        X=X, sym=sym, obs=np.asarray(a.obs_names, str), coords=coords))
        log(f"[load COMP] {sec:12s} {X.shape[0]:5d} spots {X.shape[1]:6d} symbols coords={coords is not None}")
    return out


# ---------------------------------------------------------------- load + gene space
S = load_cs() + load_comparators() + load_normal()
for s in S:
    assert abs(s["X"].data - np.round(s["X"].data)).max() < 1e-6, f"{s['section']}: non-integer counts"

union = set()
for s in S:
    union |= set(s["sym"])
inner = set(union)
for s in S:
    inner &= set(s["sym"])
inner = sorted(inner)
log(f"[join] union {len(union)} | inner-join {len(inner)}")

# gene-space accounting: which sources lack each lost gene
lost = sorted(union - set(inner))
present = {s["section"]: set(s["sym"]) for s in S}
src_of = {s["section"]: s["source"] for s in S}
rows = []
for g in lost:
    miss = [sec for sec in present if g not in present[sec]]
    msrc = sorted({src_of[sec] for sec in miss})
    rows.append(dict(gene=g, n_sections_missing=len(miss), sources_missing=";".join(msrc),
                     sections_missing=";".join(sorted(miss))))
pd.DataFrame(rows).to_csv(f"{OUT}/lost_genes.tsv", sep="\t", index=False)

acc = pd.DataFrame([dict(section=s["section"], cohort=s["cohort"], disease=s["disease"], source=s["source"],
                         n_spots=s["X"].shape[0], n_symbols=s["X"].shape[1],
                         n_symbols_in_join=len(set(s["sym"]) & set(inner)),
                         n_symbols_lost_here=len(set(s["sym"]) - set(inner)),
                         has_coords=s["coords"] is not None) for s in S])
acc["union_size"] = len(union); acc["inner_join_size"] = len(inner); acc["n_lost_total"] = len(lost)
acc.to_csv(f"{OUT}/gene_space_accounting.tsv", sep="\t", index=False)

watch = ["KCNC1", "GJA9", "FPGT", "TNNI3K", "MLIP", "PANK1", "GJB7", "NAIP", "GJA9-MYCBP", "FPGT-TNNI3K"]
wrows = [dict(gene=g, in_union=g in union, in_inner_join=g in set(inner),
              n_sections_present=sum(g in present[sec] for sec in present),
              sections_missing=";".join(sorted(sec for sec in present if g not in present[sec])))
         for g in watch]
pd.DataFrame(wrows).to_csv(f"{OUT}/watchlist_gene_presence.tsv", sep="\t", index=False)
log("[watch]\n" + pd.DataFrame(wrows).to_string(index=False))

# restrict every section to the joined space, in fixed order
inner_idx = pd.Index(inner)
for s in S:
    s["X"] = sp.csr_matrix(s["X"][:, pd.Index(s["sym"]).get_indexer(inner)])
    s["sym"] = inner_idx
    s["total"] = np.asarray(s["X"].sum(1)).ravel()

# ---------------------------------------------------------------- per-depth pipeline
lin_genes = {L: [g for g in gs if g in set(inner)] for L, gs in LIN.items()}
json.dump({L: gs for L, gs in lin_genes.items()}, open(f"{OUT}/lineage_genes_used.json", "w"), indent=1)
lin_names = list(LIN.keys())

retention, regions, spot_tabs, lin_means = [], [], [], []
for T in T_GRID:
    log(f"=== T={T} ===")
    keepmask = {}
    blocks, metas = [], []
    for s in S:
        m = s["total"] >= T
        keepmask[s["section"]] = m
        a = ad.AnnData(sp.csr_matrix(s["X"][m]).astype(np.float32))
        sc.pp.downsample_counts(a, counts_per_cell=T, random_state=RS)
        Xd = sp.csr_matrix(a.X)
        assert np.allclose(np.asarray(Xd.sum(1)).ravel(), T), f"{s['section']}: thinning did not hit T"
        blocks.append(Xd)
        metas.append(pd.DataFrame(dict(spot=s["obs"][m], section=s["section"], cohort=s["cohort"],
                                       disease=s["disease"], source=s["source"],
                                       umi_native=s["total"][m])))
        retention.append(dict(T=T, section=s["section"], cohort=s["cohort"], disease=s["disease"],
                              source=s["source"], n_spots=len(m), n_retained=int(m.sum()),
                              pct_retained=100.0 * m.mean(), median_umi_native=float(np.median(s["total"]))))
    Xall = sp.vstack(blocks).tocsr()
    md = pd.concat(metas, ignore_index=True)
    log(f"[depth] {Xall.shape[0]} spots retained at T={T}")

    # ---- lineage scores: ONE pooled z-transform fitted once over all spots at this T
    L1 = Xall.copy(); L1.data = np.log1p(L1.data)
    gidx = {g: i for i, g in enumerate(inner)}
    scores = np.zeros((Xall.shape[0], len(lin_names)))
    for j, L in enumerate(lin_names):
        cols = [gidx[g] for g in lin_genes[L]]
        Z = np.asarray(L1[:, cols].todense())
        mu, sd = Z.mean(0), Z.std(0)
        sd[sd == 0] = 1.0
        scores[:, j] = ((Z - mu) / sd).mean(1)
    dom = np.array(lin_names)[scores.argmax(1)]
    md["dominant_lineage"] = dom
    for j, L in enumerate(lin_names):
        md[f"z_{L}"] = scores[:, j]

    # ---- CS-only granuloma exclusion (comparators and normals: none)
    md["dist_gran_mm"] = np.nan
    md["excluded_granuloma_radius"] = False
    md["px_per_um"] = np.nan
    off = 0
    for s in S:
        n = int(keepmask[s["section"]].sum())
        sl = slice(off, off + n); off += n
        if s["cohort"] != "CS" or s["coords"] is None or n < 10:
            continue
        C = s["coords"][keepmask[s["section"]]]
        d, _ = cKDTree(C).query(C, k=2)
        ppu = float(np.median(d[:, 1])) / PITCH          # pixels per micron
        gm = np.isin(dom[sl], list(GRAN_LIN))
        md.loc[md.index[sl], "px_per_um"] = ppu
        if gm.sum() and np.isfinite(ppu) and ppu > 0:
            dg, _ = cKDTree(C[gm]).query(C, k=1)
            dg_mm = dg / ppu / 1000.0
            md.loc[md.index[sl], "dist_gran_mm"] = dg_mm
            md.loc[md.index[sl], "excluded_granuloma_radius"] = dg_mm < EXCL_MM

    cmdom = md["dominant_lineage"].values == "Cardiomyocyte"
    sel = cmdom & ~md["excluded_granuloma_radius"].values
    md["cm_dominant"] = cmdom
    md["selected"] = sel

    rrows = []
    for s in S:
        d = md[md.section.values == s["section"]]
        isg = d.dominant_lineage.isin(GRAN_LIN)
        rrows.append(dict(section=s["section"], cohort=s["cohort"], disease=s["disease"], source=s["source"],
                          n_retained_at_T=len(d), n_cm_dominant=int(d.cm_dominant.sum()),
                          n_granulomatous=int(isg.sum()),
                          n_excluded_radius=int((d.cm_dominant & d.excluded_granuloma_radius).sum()),
                          n_selected=int(d.selected.sum()),
                          pct_granulomatous=100.0 * isg.mean() if len(d) else np.nan,
                          median_dist_gran_mm=float(np.nanmedian(d.dist_gran_mm)) if d.dist_gran_mm.notna().any() else np.nan))
    reg = pd.DataFrame(rrows)
    reg["T"] = T
    reg["evaluable"] = reg.n_selected >= MIN_SPOTS_EVALUABLE
    regions.append(reg)
    log("[region]\n" + reg.to_string(index=False))

    # ---- pseudobulk over SELECTED spots, per section
    secs = [s["section"] for s in S]
    pb = np.zeros((len(secs), len(inner)), dtype=np.int64)
    det = {}
    for i, sec in enumerate(secs):
        m = ((md.section.values == sec) & sel)
        if m.sum() == 0:
            continue
        Xs = sp.csr_matrix(Xall[np.where(m)[0]])
        pb[i] = np.asarray(Xs.sum(0)).ravel().astype(np.int64)
        det[sec] = np.asarray((Xs > 0).sum(0)).ravel() / m.sum()
    pd.DataFrame(pb, index=secs, columns=inner).to_csv(f"{OUT}/pseudobulk_T{T}.tsv.gz", sep="\t")
    pd.DataFrame(det, index=inner).T.to_csv(f"{OUT}/detection_by_section_T{T}.tsv.gz", sep="\t")

    # group-level detection fraction over selected spots
    grp = {}
    for coh in ["CS", "COMPARATOR", "NORMAL"]:
        m = (md.cohort.values == coh) & sel
        if m.sum() == 0:
            continue
        Xs = sp.csr_matrix(Xall[np.where(m)[0]])
        grp[coh] = np.asarray((Xs > 0).sum(0)).ravel() / m.sum()
    m = sel & np.isin(md.cohort.values, ["CS", "COMPARATOR"])
    grp["DISEASED"] = np.asarray((sp.csr_matrix(Xall[np.where(m)[0]]) > 0).sum(0)).ravel() / m.sum()
    grp["ALL"] = np.asarray((sp.csr_matrix(Xall[np.where(sel)[0]]) > 0).sum(0)).ravel() / sel.sum()
    pd.DataFrame(grp, index=inner).to_csv(f"{OUT}/detection_fraction_T{T}.tsv.gz", sep="\t")

    lm = (md[sel].groupby(["section", "cohort", "disease"], observed=True)[[f"z_{L}" for L in lin_names]]
            .mean().reset_index())
    lm["T"] = T
    lin_means.append(lm)

    cols = ["spot", "section", "cohort", "disease", "source", "umi_native", "dominant_lineage",
            "dist_gran_mm", "excluded_granuloma_radius", "cm_dominant", "selected"] + [f"z_{L}" for L in lin_names]
    t = md[cols].copy(); t["T"] = T
    spot_tabs.append(t)

pd.DataFrame(retention).to_csv(f"{OUT}/retention_by_group.tsv", sep="\t", index=False)
pd.concat(regions, ignore_index=True).to_csv(f"{OUT}/region_selection.tsv", sep="\t", index=False)
pd.concat(lin_means, ignore_index=True).to_csv(f"{OUT}/section_lineage_means.tsv", sep="\t", index=False)
pd.concat(spot_tabs, ignore_index=True).to_csv(f"{OUT}/per_spot_regions.tsv.gz", sep="\t", index=False)
pd.Series(inner).to_csv(f"{OUT}/inner_join_genes_genomewide.txt", index=False, header=False)
json.dump(dict(n_sections=len(S), inner_join=len(inner), union=len(union), n_lost=len(lost),
               min_spots_evaluable=MIN_SPOTS_EVALUABLE, T_grid=T_GRID),
          open(f"{OUT}/export_summary_spatialfirst.json", "w"), indent=1)
log("[done]")
