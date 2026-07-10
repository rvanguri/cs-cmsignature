#!/usr/bin/env python3
"""Spatial visualization of key CS marker programs in the Foong CS Visium data (GEO GSE314910).

GSE314910 = 12 CS myocardial Visium samples (9 patients); CS-only, so there is no in-series CS-vs-HCM
contrast — the value here is showing WHERE the validated CS programs sit in intact tissue, especially
tertiary-lymphoid-structure (TLS)/B-cell niches. For each sample we render, per marker and per program
signature score: spatial feature maps on the tissue, plus a cross-sample expression summary.

GEO deposits Space Ranger outputs FLAT (per-sample *_filtered_feature_bc_matrix.h5 alongside a
tissue_positions CSV, scalefactors JSON, and hires PNG), not as ready Space Ranger dirs, so we
reconstruct each Visium AnnData from those flat files. Falls back to a combined .h5ad if one is present.
"""
from __future__ import annotations
import argparse, glob, json, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import log, ensure_dir  # noqa: E402
LOG = log("foong_fig")

PANEL = {
    "CS_cardiomyocyte":  ["GJB7", "TNNI3K", "MLIP", "PANK1"],
    # Granuloma-inflammation reference (matches cm_spatial_crossdisease.CS_CM_SIGS["inflammatory"];
    # 4th gene is IL7, not the deprecated CASP4). Curated, not data-derived; used as a contrast axis.
    "Inflammasome":      ["NLRC4", "IL1RAP", "BACH2", "IL7"],
    "Arrhythmogenic":    ["PERP", "CALM2", "GJA1", "SCN5A"],
    "TLS_Bcell":         ["MS4A1", "CD79A", "EBF1", "BANK1", "MZB1"],
    "TLS_chemokine":     ["CXCL13", "CCL19", "CCL21", "CR2", "LTB"],
    "Granuloma_myeloid": ["CD163", "GPNMB", "CHI3L1"],
}
# signatures to score+map spatially (the TLS score is the headline for the letter)
SIGNATURES = {
    "TLS_score": PANEL["TLS_Bcell"] + PANEL["TLS_chemokine"],
    "Inflammasome_score": PANEL["Inflammasome"],
    "CScardiomyocyte_score": PANEL["CS_cardiomyocyte"],
}


def _open_maybe_gz(path, text=True):
    import gzip
    if path.endswith(".gz"):
        return gzip.open(path, "rt" if text else "rb")
    return open(path, "r" if text else "rb")


def _read_positions(files):
    import pandas as pd
    cand = [f for f in files if re.search(r"tissue_positions.*\.csv(\.gz)?$", f)] or \
           [f for f in files if re.search(r"tissue_positions.*\.parquet$", f)]
    if not cand:
        return None
    p = cand[0]
    if p.endswith(".parquet"):
        df = pd.read_parquet(p)
    else:
        with _open_maybe_gz(p) as fh:      # pandas reads .gz natively, but sniff the header first
            head = fh.readline()
        if "barcode" in head or "in_tissue" in head:
            df = pd.read_csv(p)            # has header (SpaceRanger >=2.0); read_csv handles .gz
        else:                             # older tissue_positions_list.csv: no header
            df = pd.read_csv(p, header=None, names=["barcode", "in_tissue", "array_row",
                             "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"])
    df.columns = [c.lower() for c in df.columns]
    return df.set_index("barcode")


def _load_scalefactors(files):
    sf = [f for f in files if re.search(r"scalefactors_json.*\.json(\.gz)?$", f)]
    if not sf:
        return None
    try:
        with _open_maybe_gz(sf[0]) as fh:
            return json.load(fh)
    except Exception as e:  # noqa: BLE001
        LOG.warning("scalefactors load failed (%s)", e); return None


def _load_hires(files):
    img = [f for f in files if re.search(r"tissue_hires_image.*\.png(\.gz)?$", f)]
    if not img:
        return None
    try:
        import io, matplotlib.image as mpimg
        with _open_maybe_gz(img[0], text=False) as fh:
            return mpimg.imread(io.BytesIO(fh.read()), format="png")
    except Exception as e:  # noqa: BLE001
        LOG.warning("hires image load failed (%s)", e); return None


def plot_pair(a, sample, color, title, out_png):
    """Side-by-side: left = clean H&E (no dots), right = same H&E with the expression overlay.

    Both panels are framed identically (the H&E is cropped to the overlay's axis limits). If no
    tissue image is available we fall back to a single spot-scatter panel.
    """
    import matplotlib.pyplot as plt
    import scanpy as sc
    img = None
    if "spatial" in a.uns:
        img = a.uns["spatial"].get(sample, {}).get("images", {}).get("hires")
    if img is not None:
        fig, (ax_he, ax_ov) = plt.subplots(1, 2, figsize=(12, 6))
        sc.pl.spatial(a, color=color, ax=ax_ov, title=title, show=False)   # H&E background + dots
        ax_he.imshow(img)
        ax_he.set_xlim(ax_ov.get_xlim()); ax_he.set_ylim(ax_ov.get_ylim())  # match framing
        ax_he.set_title(f"{sample} H&E"); ax_he.axis("off")
        fig.tight_layout(); fig.savefig(out_png, dpi=150, bbox_inches="tight"); plt.close(fig)
    else:
        spot = float(a.uns.get("_sf", {}).get("spot_diameter_fullres", 100.0))
        sc.pl.spatial(a, color=color, img_key=None, spot_size=spot, title=title, show=False)
        plt.savefig(out_png, dpi=150, bbox_inches="tight"); plt.close()


def load_geo_visium(root):
    """Reconstruct (sample, AnnData) from GEO's flat Space Ranger files, with spatial coords.

    GEO (GSE314910) names files GSM<id>_<descriptor>_<SAMPLE>.<ext>[.gz], e.g.
    GSM9416129_filtered_feature_bc_matrix_CS_101-1.h5 and GSM9416129_tissue_positions_CS_101-1.csv.gz,
    so we group siblings by the shared GSM<id> prefix (the sample name trails the descriptor).
    """
    import scanpy as sc
    h5ads = glob.glob(os.path.join(root, "**", "*.h5ad"), recursive=True)
    if h5ads:
        import anndata as ad
        return [(os.path.basename(f).replace(".h5ad", ""), ad.read_h5ad(f)) for f in h5ads]

    mats = sorted(glob.glob(os.path.join(root, "**", "*filtered_feature_bc_matrix*.h5"), recursive=True))
    allfiles = glob.glob(os.path.join(root, "**", "*"), recursive=True)
    out = []
    for m in mats:
        d = os.path.dirname(m); base = os.path.basename(m)
        gsm_match = re.match(r"(GSM\d+)_", base)
        gsm = gsm_match.group(1) if gsm_match else base.split("_")[0]
        sample = re.sub(r"\.h5$", "", re.sub(r"^GSM\d+_filtered_feature_bc_matrix[_-]?", "", base)) or gsm
        sibs = [f for f in allfiles if os.path.dirname(f) == d and os.path.basename(f).startswith(gsm + "_")]
        try:
            a = sc.read_10x_h5(m); a.var_names_make_unique()
        except Exception as e:  # noqa: BLE001
            LOG.warning("read_10x_h5 failed for %s (%s)", m, e); continue
        pos = _read_positions(sibs)
        if pos is not None:
            common = a.obs_names.intersection(pos.index)
            a = a[common].copy(); pos = pos.reindex(a.obs_names)
            if "in_tissue" in pos:
                a = a[pos["in_tissue"].astype(int) == 1].copy(); pos = pos.reindex(a.obs_names)
            a.obsm["spatial"] = pos[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy()
            scal = _load_scalefactors(sibs)
            if scal is not None:
                a.uns["_sf"] = scal            # keep for spot_size fallback
                hires = _load_hires(sibs)
                if hires is not None:
                    a.uns["spatial"] = {sample: {"images": {"hires": hires}, "scalefactors": scal}}
        else:
            LOG.warning("%s: no tissue_positions found -> maps will use spot scatter only", sample)
        a.obs["sample"] = sample
        LOG.info("  reconstructed %s: %d spots x %d genes (spatial=%s, image=%s)",
                 sample, a.n_obs, a.n_vars, "spatial" in a.obsm, "spatial" in a.uns)
        out.append((sample, a))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True, help="Foong download dir OR a single .h5ad")
    ap.add_argument("--out", required=True)
    ap.add_argument("--sample_map", default=None, help="optional TSV sample<TAB>label (e.g. add HCM samples)")
    args = ap.parse_args()
    ensure_dir(args.out)

    import numpy as np
    import pandas as pd
    import scanpy as sc
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    samples = load_geo_visium(args.data) if os.path.isdir(args.data) else \
        [(os.path.basename(args.data).replace(".h5ad", ""), __import__("anndata").read_h5ad(args.data))]
    if not samples:
        sys.exit("no readable Visium/h5ad data under %s" % args.data)
    LOG.info("loaded %d sample(s): %s", len(samples), [s for s, _ in samples])

    sumrows = []
    for name, a in samples:
        a.var_names = a.var_names.astype(str).str.upper(); a.var_names_make_unique()
        if float(np.asarray(a.X.max())) > 50:           # raw counts -> normalize
            sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
        # program signature scores
        for sig, genes in SIGNATURES.items():
            gs = [g for g in genes if g in a.var_names]
            if gs:
                sc.tl.score_genes(a, gs, score_name=sig)

        # ---- spatial maps: each is H&E (clean) | expression overlay ----
        sdir = os.path.join(args.out, "spatial", name); ensure_dir(sdir)
        if "spatial" in a.obsm:
            for sig in SIGNATURES:                        # signature scores (headline: TLS_score)
                if sig in a.obs:
                    try:
                        plot_pair(a, name, sig, f"{name}:{sig}", os.path.join(sdir, f"SIG_{sig}.png"))
                    except Exception as e:  # noqa: BLE001
                        LOG.warning("spatial sig %s/%s failed (%s)", name, sig, e)
            for grp, genes in PANEL.items():              # individual markers
                for g in genes:
                    if g not in a.var_names:
                        continue
                    try:
                        plot_pair(a, name, g, f"{name}:{g}", os.path.join(sdir, f"{grp}_{g}.png"))
                    except Exception as e:  # noqa: BLE001
                        LOG.warning("spatial %s/%s failed (%s)", name, g, e)
        else:
            LOG.info("%s: no spatial coords -> skipping tissue maps", name)

        # ---- per-sample mean expression of the panel ----
        flat = [g for gs in PANEL.values() for g in gs if g in a.var_names]
        mean = sc.get.obs_df(a, keys=flat).mean()
        mean["sample"] = name
        for sig in SIGNATURES:
            if sig in a.obs:
                mean[sig] = float(a.obs[sig].mean())
        sumrows.append(mean)

    # ---- cross-sample summary table + heatmap ----
    summ = pd.DataFrame(sumrows).set_index("sample")
    summ.to_csv(os.path.join(args.out, "foong_marker_summary_by_sample.tsv"), sep="\t")
    try:
        import matplotlib.pyplot as plt
        genecols = [c for c in summ.columns if c not in SIGNATURES]
        fig, ax = plt.subplots(figsize=(max(6, len(genecols) * 0.35), max(3, len(summ) * 0.35)))
        im = ax.imshow(summ[genecols].to_numpy(), aspect="auto", cmap="viridis")
        ax.set_xticks(range(len(genecols))); ax.set_xticklabels(genecols, rotation=90, fontsize=6)
        ax.set_yticks(range(len(summ))); ax.set_yticklabels(summ.index, fontsize=7)
        fig.colorbar(im, ax=ax, label="mean log-norm expr")
        fig.savefig(os.path.join(args.out, "foong_marker_heatmap.png"), dpi=200, bbox_inches="tight"); plt.close()
    except Exception as e:  # noqa: BLE001
        LOG.warning("heatmap failed (%s)", e)
    LOG.info("figures written -> %s", args.out)


if __name__ == "__main__":
    main()
