#!/usr/bin/env python3
"""Within-Foong regional analysis of the CS-CM signature (no cross-dataset batch confound).

Entirely inside the CS Visium data (Foong GSE314910). Assigns each spot a histologic-like ZONE from
expression: preserved (cardiomyocyte-dominant), granulomatous (immune/myeloid-dominant), or fibrotic
(fibroblast/ECM-dominant): then asks WHERE the two signatures live:

  - CM_intrinsic (GJB7/TNNI3K/MLIP/PANK1) and inflammatory (NLRC4/IL1RAP/BACH2/IL7) score per zone
    (violin + Kruskal-Wallis + preserved-vs-lesional Mann-Whitney).
  - Distance-to-lesion curve: for preserved spots, is the CM_intrinsic score flat with distance from
    the nearest lesional (granulomatous/fibrotic) spot (cell-autonomous) or does it rise toward
    lesions (paracrine)?  Reported as a within-patient linear mixed model (spots nested in patients) + per-patient slope signs.
  - Per-sample spatial zone maps and CM_intrinsic maps (H&E-paired).

Reuses the Visium loader + lineage scoring from the cross-disease script.
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import log, ensure_dir  # noqa: E402
from foong_spatial_figures import load_geo_visium, plot_pair  # noqa: E402
from cm_spatial_crossdisease import ensure_symbols, LINEAGES, CS_CM_SIGS  # noqa: E402
LOG = log("foong_region")

# zone <- dominant lineage program
ZONE_OF = {"Cardiomyocyte": "preserved", "Myeloid": "granulomatous", "Lymphoid": "granulomatous",
           "Fibroblast": "fibrotic", "Endothelial": "other"}
LESIONAL = {"granulomatous", "fibrotic"}


def assign_zones(a, sc):
    lin = {}
    for L, gs in LINEAGES.items():
        gg = [g for g in gs if g in a.var_names]
        if gg:
            sc.tl.score_genes(a, gg, score_name=f"_l_{L}", use_raw=False)
            lin[L] = a.obs[f"_l_{L}"].values
    if "Cardiomyocyte" not in lin:
        return None
    dom = pd.DataFrame(lin).idxmax(axis=1).values
    return np.array([ZONE_OF.get(d, "other") for d in dom])


def nearest_lesion_dist(coords, is_lesion):
    """Euclidean distance from each spot to the nearest lesional spot (same sample)."""
    from scipy.spatial import cKDTree
    if is_lesion.sum() == 0 or (~is_lesion).sum() == 0:
        return np.full(len(coords), np.nan)
    tree = cKDTree(coords[is_lesion])
    d, _ = tree.query(coords, k=1)
    return d


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=None, help="Foong Visium dir (defaults to $RAW/foong_spatial)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--n_dist_bins", type=int, default=8)
    args = ap.parse_args()
    data = args.data or os.path.join(os.environ.get("RAW", "."), "foong_spatial")
    ensure_dir(args.out)

    import scanpy as sc
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.stats import kruskal, mannwhitneyu, spearmanr

    samples = load_geo_visium(data)
    if not samples:
        sys.exit("no Visium samples under %s" % data)
    LOG.info("loaded %d CS samples", len(samples))

    spot_rows, dist_rows = [], []
    for name, a in samples:
        a = ensure_symbols(a)
        a.var_names = a.var_names.astype(str).str.upper(); a.var_names_make_unique()
        if float(np.asarray(a.X.max())) > 50:
            sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
        zones = assign_zones(a, sc)
        if zones is None:
            LOG.warning("%s: no CM markers -> skip", name); continue
        a.obs["zone"] = zones
        for sig, genes in CS_CM_SIGS.items():
            gg = [g for g in genes if g in a.var_names]
            if gg:
                sc.tl.score_genes(a, gg, score_name=sig, use_raw=False)
        # per-spot record
        for i in range(a.n_obs):
            r = {"sample": name, "zone": zones[i]}
            for sig in CS_CM_SIGS:
                if sig in a.obs:
                    r[sig] = float(a.obs[sig].values[i])
            spot_rows.append(r)
        # distance-to-lesion for preserved spots
        if "spatial" in a.obsm and "CM_intrinsic" in a.obs:
            coords = np.asarray(a.obsm["spatial"], float)
            is_les = np.isin(zones, list(LESIONAL))
            d = nearest_lesion_dist(coords, is_les)
            pres = zones == "preserved"
            for i in np.where(pres & np.isfinite(d))[0]:
                dist_rows.append({"sample": name, "dist": float(d[i]),
                                  "CM_intrinsic": float(a.obs["CM_intrinsic"].values[i])})
        # zone + CM_intrinsic spatial maps
        if "spatial" in a.obsm:
            try:
                a.obs["zone_code"] = pd.Categorical(zones).codes.astype(float)
                plot_pair(a, name, "CM_intrinsic", f"{name}: CM-intrinsic",
                          os.path.join(args.out, f"spatial_{name}_CMintrinsic.png"))
            except Exception as e:  # noqa: BLE001
                LOG.warning("spatial %s failed (%s)", name, e)
        n_zone = pd.Series(zones).value_counts().to_dict()
        LOG.info("%s: zones %s", name, n_zone)

    df = pd.DataFrame(spot_rows)
    df.to_csv(os.path.join(args.out, "foong_regional_perspot.tsv.gz"), sep="\t", index=False)

    # ---- per-zone summary: MEAN and median per signature, each zone broken out separately ----
    sig_cols = [s for s in CS_CM_SIGS if s in df.columns]
    zsummary = df.groupby("zone")[sig_cols].agg(["mean", "median", "size"])
    zsummary.to_csv(os.path.join(args.out, "foong_zone_means.tsv"), sep="\t")
    LOG.info("per-zone MEAN of each signature:\n%s",
             df.groupby("zone")[sig_cols].mean().round(4).to_string())

    # ---- signature by zone: violin + stats ----
    zorder = [z for z in ["preserved", "granulomatous", "fibrotic", "other"] if z in df.zone.unique()]
    stats = []
    for sig in [s for s in CS_CM_SIGS if s in df.columns]:
        groups = [df.loc[df.zone == z, sig].dropna().values for z in zorder]
        fig, ax = plt.subplots(figsize=(1.6 * len(zorder) + 2, 5))
        ax.violinplot(groups, showmeans=True, showextrema=False)
        ax.boxplot(groups, widths=0.15, showfliers=False)
        ax.set_xticks(range(1, len(zorder) + 1)); ax.set_xticklabels(zorder, rotation=20)
        ax.set_ylabel(f"{sig} score"); ax.set_title(f"{sig} by CS tissue zone (within-Foong)")
        fig.tight_layout(); fig.savefig(os.path.join(args.out, f"zone_{sig}.png"), dpi=200,
                                        bbox_inches="tight"); plt.close(fig)
        try:
            _, p_kw = kruskal(*[g for g in groups if len(g)])
        except Exception:  # noqa: BLE001
            p_kw = np.nan
        pres = df.loc[df.zone == "preserved", sig].dropna()
        gran = df.loc[df.zone == "granulomatous", sig].dropna()   # break out separately (not pooled)
        fibr = df.loc[df.zone == "fibrotic", sig].dropna()
        def _mwu(a, b):
            try:
                return mannwhitneyu(a, b, alternative="two-sided")[1]
            except Exception:  # noqa: BLE001
                return np.nan
        stats.append({
            "signature": sig,
            "mean_preserved": float(pres.mean()) if len(pres) else np.nan,
            "mean_granulomatous": float(gran.mean()) if len(gran) else np.nan,
            "mean_fibrotic": float(fibr.mean()) if len(fibr) else np.nan,
            "kruskal_p_across_zones": p_kw,
            "mwu_p_granuloma_vs_preserved": _mwu(gran, pres),
        })
    pd.DataFrame(stats).to_csv(os.path.join(args.out, "foong_zone_stats.tsv"), sep="\t", index=False)

    # ---- distance-to-lesion: WITHIN-PATIENT mixed-effects model ----
    # A pooled Spearman across spots is pseudoreplicated (thousands of correlated spots per patient)
    # and vulnerable to Simpson's paradox. The primary statistic is a linear mixed model with spots
    # nested in patients (random intercept); we also report the per-patient slope signs. Distance is
    # z-scored WITHIN each sample so tissues of different physical scale are comparable.
    dd = pd.DataFrame(dist_rows)
    if len(dd) > 50 and dd["sample"].nunique() >= 3:
        dd["patient"] = dd["sample"].astype(str).str.replace(r"-\d+$", "", regex=True)
        dd["dist_z"] = dd.groupby("sample")["dist"].transform(lambda x: (x - x.mean()) / (x.std() + 1e-9))
        dd = dd.replace([np.inf, -np.inf], np.nan).dropna(subset=["dist_z", "CM_intrinsic"])
        dd.to_csv(os.path.join(args.out, "foong_distance_to_lesion.tsv.gz"), sep="\t", index=False)
        n_spots, n_pat = len(dd), dd["patient"].nunique()
        row = {"n_spots": n_spots, "n_patients": n_pat}
        try:
            import statsmodels.formula.api as smf
            m = smf.mixedlm("CM_intrinsic ~ dist_z", dd, groups=dd["patient"].values).fit(reml=True)
            beta = float(m.params["dist_z"]); ci = m.conf_int().loc["dist_z"].values
            row.update(mixedlm_beta_per_SD=beta, ci_low=float(ci[0]), ci_high=float(ci[1]),
                       mixedlm_p=float(m.pvalues["dist_z"]))
        except Exception as e:  # noqa: BLE001
            LOG.warning("mixedlm failed (%s): falling back to per-patient slopes only", e)
            row.update(mixedlm_beta_per_SD=np.nan, ci_low=np.nan, ci_high=np.nan, mixedlm_p=np.nan)
        # per-patient independent slopes (robustness: is the sign consistent?)
        slopes = [np.polyfit(g["dist_z"], g["CM_intrinsic"], 1)[0]
                  for _, g in dd.groupby("patient") if g["dist_z"].nunique() > 2]
        row["patients_positive_slope"] = int(np.sum(np.array(slopes) > 0))
        row["patients_tested"] = len(slopes)
        rho, p_sp = spearmanr(dd["dist_z"], dd["CM_intrinsic"])   # kept only as a secondary descriptor
        row.update(pooled_spearman_rho=float(rho), pooled_spearman_p=float(p_sp))
        row["interpretation"] = ("beta>=0 / positive within-patient slopes => NOT enriched toward "
                                 "lesions => inconsistent with paracrine induction")
        pd.DataFrame([row]).to_csv(os.path.join(args.out, "foong_distance_stats.tsv"), sep="\t", index=False)
        # figure: decile means + spot-level fit, annotated with the mixed-model result
        dd["bin"] = pd.qcut(dd["dist_z"].rank(method="first"), args.n_dist_bins, labels=False)
        curve = dd.groupby("bin").agg(dist=("dist_z", "mean"), score=("CM_intrinsic", "mean"),
                                      sem=("CM_intrinsic", lambda x: x.std() / max(1, np.sqrt(len(x))))).reset_index()
        fig, ax = plt.subplots(figsize=(6, 5))
        ax.errorbar(curve["dist"], curve["score"], yerr=curve["sem"], marker="o", ls="none",
                    capsize=3, color="#333", label=f"decile mean ± SEM ({args.n_dist_bins} bins)")
        b = np.polyfit(dd["dist_z"], dd["CM_intrinsic"], 1)
        xs = np.linspace(dd["dist_z"].min(), dd["dist_z"].max(), 40)
        ax.plot(xs, b[0] * xs + b[1], "r-", lw=1.5, label="spot-level linear fit")
        ax.set_xlabel("distance to nearest lesion (within-sample z)")
        ax.set_ylabel("CM-intrinsic score (preserved spots)")
        # report a P THRESHOLD, not the raw exponent: with ~10^4 spots even a near-zero slope is
        # astronomically 'significant', so an extreme exponent misleads. The claim is direction + 8/8.
        _p = row.get("mixedlm_p", float("nan"))
        _ptxt = "P<0.001" if (_p == _p and _p < 1e-3) else f"P={_p:.2g}"
        ax.set_title(f"shallow positive slope β={row.get('mixedlm_beta_per_SD', float('nan')):+.3f}/SD ({_ptxt})\n"
                     f"positive in {row['patients_positive_slope']}/{row['patients_tested']} patients "
                     f": direction, not magnitude")
        ax.legend(fontsize=7, frameon=False)
        fig.tight_layout(); fig.savefig(os.path.join(args.out, "distance_to_lesion.png"), dpi=200,
                                        bbox_inches="tight"); plt.close(fig)
        LOG.info("distance mixed model: β=%.4f p=%.2e  positive in %d/%d patients  (n=%d spots, %d patients)",
                 row.get("mixedlm_beta_per_SD", float("nan")), row.get("mixedlm_p", float("nan")),
                 row["patients_positive_slope"], row["patients_tested"], n_spots, n_pat)
    else:
        LOG.warning("too few preserved spots/patients for the lesion mixed model")

    LOG.info("zone counts:\n%s", df.groupby("zone").size().to_string())
    LOG.info("done -> %s", args.out)


if __name__ == "__main__":
    main()
