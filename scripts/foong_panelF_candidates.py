#!/usr/bin/env python3
"""Panel-F candidates: for EACH Foong CS Visium sample, a 3-panel strip:
   [ clean H&E | CM-intrinsic score | inferred granuloma (immune-dominant zone) ]
so the best representative sample for the figure can be chosen by eye.

Granuloma annotation is COMPUTATIONAL (immune/myeloid-dominant spots), not pathologist-drawn; labeled
as such. Runs on the cluster (needs the Visium data, default $RAW/foong_spatial).
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import log, ensure_dir  # noqa: E402
from foong_spatial_figures import load_geo_visium  # noqa: E402
from cm_spatial_crossdisease import ensure_symbols, CS_CM_SIGS, LINEAGES  # noqa: E402
from foong_regional import assign_zones  # noqa: E402
LOG = log("panelF")

ZONE_ORDER = ["preserved", "granulomatous", "fibrotic", "other"]
ZONE_COLORS = ["#d9d9d9", "#d62728", "#ff7f0e", "#f2f2f2"]   # granuloma = red


def _score_map(ax, xh, yh, vals, img, cmap, mode, label):
    """One spatial score panel: H&E background + score-colored spots (alpha).
    mode: 'sequential' (p2-p98), 'high' (p55-p98: emphasize the high end, quiet tissue -> dark),
          'diverging' (symmetric around 0)."""
    import matplotlib.pyplot as plt
    vals = np.asarray(vals, float)
    if img is not None:
        ax.imshow(img)
    if mode == "diverging":
        vmax = float(np.nanpercentile(np.abs(vals), 98)) or 1e-6
        kw = dict(cmap=cmap, vmin=-vmax, vmax=vmax)
    elif mode == "high":
        kw = dict(cmap=cmap, vmin=float(np.nanpercentile(vals, 55)),
                  vmax=float(np.nanpercentile(vals, 98)))
    else:
        kw = dict(cmap=cmap, vmin=float(np.nanpercentile(vals, 2)),
                  vmax=float(np.nanpercentile(vals, 98)))
    sm = ax.scatter(xh, yh, c=vals, s=7, alpha=0.8, linewidths=0, **kw)
    plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.02, label=label)
    ax.axis("off")


def single_panel(a, name, out_png):
    """Figure-ready single panel: H&E + CM-intrinsic score spots + RED contour of the inferred
    (immune-dominant) granuloma. One panel -> respects the 6-panel Research Letter limit."""
    import matplotlib.pyplot as plt
    from scipy.ndimage import gaussian_filter
    has_img = "spatial" in a.uns
    hscf = float(a.uns.get("_sf", {}).get("tissue_hires_scalef", 1.0)) if has_img else 1.0
    coords = np.asarray(a.obsm["spatial"], float)
    xh, yh = coords[:, 0] * hscf, coords[:, 1] * hscf
    fig, ax = plt.subplots(figsize=(6.5, 6))
    if has_img:
        img = list(a.uns["spatial"].values())[0]["images"].get("hires")
        ax.imshow(img); Wd, Ht = img.shape[1], img.shape[0]
    else:
        Wd, Ht = xh.max() * 1.05, yh.max() * 1.05
        ax.set_xlim(0, Wd); ax.set_ylim(Ht, 0)
    sctt = ax.scatter(xh, yh, c=a.obs["CM_intrinsic"].values, cmap="viridis", s=6, alpha=0.75, linewidths=0)
    plt.colorbar(sctt, ax=ax, fraction=0.046, pad=0.02, label="CM-intrinsic score")
    is_g = a.obs["zone"].astype(str).values == "granulomatous"
    if int(is_g.sum()) >= 5:                                   # red outline = granuloma-spot density contour
        Hh, xe, ye = np.histogram2d(xh[is_g], yh[is_g], bins=70, range=[[0, Wd], [0, Ht]])
        Hh = gaussian_filter(Hh.T, sigma=1.4)
        if Hh.max() > 0:
            ax.contour(0.5 * (xe[:-1] + xe[1:]), 0.5 * (ye[:-1] + ye[1:]), Hh,
                       levels=[Hh.max() * 0.25], colors="red", linewidths=1.6)
    ax.set_title(f"{name}: CM-intrinsic score in situ\n(red outline = inferred granuloma, immune-dominant)", fontsize=9)
    ax.axis("off")
    fig.tight_layout(); fig.savefig(out_png, dpi=220, bbox_inches="tight"); plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=None, help="Foong Visium dir (default $RAW/foong_spatial)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    data = args.data or os.path.join(os.environ.get("RAW", "."), "foong_spatial")
    ensure_dir(args.out)

    import scanpy as sc
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    samples = load_geo_visium(data)
    if not samples:
        sys.exit("no Visium samples under %s" % data)
    LOG.info("loaded %d samples", len(samples))
    summary = []

    for name, a in samples:
        a = ensure_symbols(a)
        a.var_names = a.var_names.astype(str).str.upper(); a.var_names_make_unique()
        if float(np.asarray(a.X.max())) > 50:
            sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
        zones = assign_zones(a, sc)
        if zones is None:
            LOG.warning("%s: no CM markers -> skip", name); continue
        a.obs["zone"] = pd.Categorical(zones, categories=ZONE_ORDER)
        gg = [g for g in CS_CM_SIGS["CM_intrinsic"] if g in a.var_names]
        sc.tl.score_genes(a, gg, score_name="CM_intrinsic", use_raw=False)
        # continuous granuloma score = immune (myeloid + lymphoid) program
        imm = [g for L in ("Myeloid", "Lymphoid") for g in LINEAGES[L] if g in a.var_names]
        sc.tl.score_genes(a, imm, score_name="granuloma_score", use_raw=False)

        has_img = "spatial" in a.uns
        img = list(a.uns["spatial"].values())[0]["images"].get("hires") if has_img else None
        hscf = float(a.uns.get("_sf", {}).get("tissue_hires_scalef", 1.0)) if has_img else 1.0
        coords = np.asarray(a.obsm["spatial"], float)
        xh, yh = coords[:, 0] * hscf, coords[:, 1] * hscf

        # Panel F: H&E | granuloma (immune) score | CM-intrinsic score.
        # Colormaps per spatial-transcriptomics convention: sequential 'magma' for the one-directional
        # granuloma score; DIVERGING 'coolwarm' (centered at 0) for the CM-intrinsic score, which is a
        # score_genes value with a meaningful zero: easier to read than viridis, with alpha over H&E.
        fig, ax = plt.subplots(1, 3, figsize=(16, 5.2))
        if img is not None:
            ax[0].imshow(img)
        ax[0].set_title(f"{name}: H&E", fontsize=10); ax[0].axis("off")
        _score_map(ax[1], xh, yh, a.obs["granuloma_score"].values, img, "magma", "sequential",
                   "immune-cell density")
        ax[1].set_title(f"{name}: immune-cell density (granuloma)", fontsize=10)
        # CM-intrinsic: sequential 'plasma' clipped to the high end (quiet tissue -> dark) reads far
        # more clearly than a diverging map whose white centre blends into the H&E.
        _score_map(ax[2], xh, yh, a.obs["CM_intrinsic"].values, img, "plasma", "high",
                   "CM-intrinsic score")
        ax[2].set_title(f"{name}: CM-intrinsic score", fontsize=10)
        if img is not None:                                    # frame H&E to match the score panels
            ax[0].set_xlim(ax[2].get_xlim()); ax[0].set_ylim(ax[2].get_ylim())

        comp = pd.Series(zones).value_counts()
        frac = {z: round(100 * comp.get(z, 0) / len(zones), 1) for z in ZONE_ORDER}
        fig.suptitle(f"{name} :  zones: preserved {frac['preserved']}%, "
                     f"granulomatous {frac['granulomatous']}%, fibrotic {frac['fibrotic']}%  "
                     f"(n={len(zones)} spots)", fontsize=9, y=1.02)
        fig.tight_layout()
        fig.savefig(os.path.join(args.out, f"panelF_{name}.png"), dpi=160, bbox_inches="tight")
        plt.close(fig)
        # figure-ready single panel (H&E + score + red granuloma outline)
        try:
            single_panel(a, name, os.path.join(args.out, f"panelF_single_{name}.png"))
        except Exception as e:  # noqa: BLE001
            LOG.warning("%s single panel failed (%s)", name, e)
        summary.append({"sample": name, "n_spots": len(zones), **frac})
        LOG.info("%s zones %% %s", name, frac)

    pd.DataFrame(summary).to_csv(os.path.join(args.out, "panelF_zone_summary.tsv"), sep="\t", index=False)
    LOG.info("done -> %s  (one panelF_<sample>.png per sample + panelF_zone_summary.tsv)", args.out)
    LOG.info("pick a sample with substantial preserved AND granulomatous zones and clean tissue.")


if __name__ == "__main__":
    main()
