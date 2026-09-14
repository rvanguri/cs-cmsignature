#!/usr/bin/env python3
"""Gene-set robustness of the two spatial CS conclusions (cell-autonomy + stage-independence).

Re-scores the cardiomyocyte program on DEFENSIBLE gene sets instead of the hand-curated four genes,
then re-runs the IDENTICAL published analyses per gene set. The analysis logic is imported verbatim
from foong_regional.py / cm_spatial_crossdisease.py; the ONLY thing that varies across runs is which
gene list is scored into the "CM program" column.

Gene sets (metadata/geneset_*.txt, one symbol per line):
  original4    GJB7,TNNI3K,MLIP,PANK1            -- the curated panel (baseline)
  objective53  53 genes                          -- pre-specified objective rule
  tnni3k       TNNI3K only                       -- single standout gene (raw expr, not score_genes)
  mlip_pank1   MLIP,PANK1                         -- the only batch-clean (Liu CS>ICM) pair
  negctrl53    53 random CM-expressed genes       -- NEGATIVE CONTROL (must NOT reproduce the pattern)

Analyses (each emits a tidy table with a `geneset` column):
  A  distance-to-lesion mixed model   (cell-autonomy / anti-paracrine)   -> Foong Visium
  B  zone gradient                    (stage-independence)               -> Foong Visium
  C  snRNA zone flatness              (stage-independence cross-check)    -> CM pseudobulk (no raw data)
  D  cross-disease ordering           (CS vs normal vs mimics)           -> cross-disease manifest

Single-gene sets: score_genes is undefined for n=1, so the gene's log-normalized expression is used
directly (documented caveat: methodological difference from the multi-gene module score).
"""
from __future__ import annotations
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import log, ensure_dir  # noqa: E402
from foong_spatial_figures import load_geo_visium  # noqa: E402
from cm_spatial_crossdisease import ensure_symbols, LINEAGES  # noqa: E402
from foong_regional import assign_zones, nearest_lesion_dist, LESIONAL  # noqa: E402

LOG = log("geneset_robust")

ZONE_ORDER = ["preserved", "granulomatous", "fibrotic", "other"]


# ---------------------------------------------------------------------------
# gene-set IO + scoring
# ---------------------------------------------------------------------------
def read_geneset(path):
    """One symbol per line; blank lines / comments ignored. Returns (label, [genes])."""
    label = os.path.basename(path)
    label = label[len("geneset_"):] if label.startswith("geneset_") else label
    label = os.path.splitext(label)[0]
    genes = []
    with open(path) as fh:
        for line in fh:
            g = line.strip().upper()
            if g and not g.startswith("#"):
                genes.append(g)
    return label, genes


def score_set(a, genes, sc, score_name):
    """Score a gene set into a.obs[score_name].

    n>=2 -> scanpy score_genes (the published module-score routine, use_raw=False, deterministic).
    n==1 -> the single gene's log-normalized expression directly (score_genes is undefined for n=1).
    Returns the list of genes actually present (so callers can report detection).
    """
    present = [g for g in genes if g in a.var_names]
    if not present:
        a.obs[score_name] = np.nan
        return present
    if len(present) == 1:
        col = sc.get.obs_df(a, keys=present)[present[0]]
        a.obs[score_name] = col.to_numpy(dtype=float)
    else:
        sc.tl.score_genes(a, present, score_name=score_name, use_raw=False, random_state=0)
    return present


# ---------------------------------------------------------------------------
# A + B : within-Foong distance-to-lesion + zone gradient (identical to foong_regional.py)
# ---------------------------------------------------------------------------
def run_within_cs(foong_dir, genesets, out):
    import scanpy as sc
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.stats import kruskal, mannwhitneyu, wilcoxon

    samples = load_geo_visium(foong_dir)
    if not samples:
        sys.exit(f"no Visium samples under {foong_dir}")
    LOG.info("loaded %d CS samples for within-CS analysis", len(samples))

    labels = [lab for lab, _ in genesets]
    # per-spot records: one row per spot per sample carries zone, distance, and every gene-set score
    spot_rows = []
    # distance records: preserved spots only, one row per spot carries every gene-set score
    dist_rows = []

    for name, a in samples:
        a = ensure_symbols(a)
        a.var_names = a.var_names.astype(str).str.upper()
        a.var_names_make_unique()
        if float(np.asarray(a.X.max())) > 50:
            sc.pp.normalize_total(a, target_sum=1e4)
            sc.pp.log1p(a)
        zones = assign_zones(a, sc)
        if zones is None:
            LOG.warning("%s: no CM markers -> skip", name)
            continue
        a.obs["zone"] = zones
        detection = {}
        for lab, genes in genesets:
            present = score_set(a, genes, sc, f"score_{lab}")
            detection[lab] = f"{len(present)}/{len(genes)}"
        LOG.info("%s: n_spots=%d zones=%s detection=%s", name, a.n_obs,
                 pd.Series(zones).value_counts().to_dict(), detection)

        for i in range(a.n_obs):
            r = {"sample": name, "zone": zones[i]}
            for lab in labels:
                r[lab] = float(a.obs[f"score_{lab}"].values[i])
            spot_rows.append(r)

        # distance-to-lesion for preserved spots (same coords/lesion definition as foong_regional)
        if "spatial" in a.obsm:
            coords = np.asarray(a.obsm["spatial"], float)
            is_les = np.isin(zones, list(LESIONAL))
            d = nearest_lesion_dist(coords, is_les)
            pres = zones == "preserved"
            for i in np.where(pres & np.isfinite(d))[0]:
                r = {"sample": name, "dist": float(d[i])}
                for lab in labels:
                    r[lab] = float(a.obs[f"score_{lab}"].values[i])
                dist_rows.append(r)

    spot_df = pd.DataFrame(spot_rows)
    spot_df.to_csv(os.path.join(out, "within_cs_perspot.tsv.gz"), sep="\t", index=False)
    dist_df = pd.DataFrame(dist_rows)

    # ---- (A) distance-to-lesion mixed model, per gene set ----
    import statsmodels.formula.api as smf
    dist_df["patient"] = dist_df["sample"].astype(str).str.replace(r"-\d+$", "", regex=True)
    dist_df["dist_z"] = dist_df.groupby("sample")["dist"].transform(
        lambda x: (x - x.mean()) / (x.std() + 1e-9))
    dist_df = dist_df.replace([np.inf, -np.inf], np.nan)

    a_rows = []
    for lab in labels:
        dd = dist_df.dropna(subset=["dist_z", lab]).copy()
        if len(dd) < 50 or dd["patient"].nunique() < 3:
            a_rows.append({"geneset": lab, "note": "too few preserved spots/patients"})
            continue
        rec = {"geneset": lab, "n_spots": len(dd), "n_patients": dd["patient"].nunique()}
        try:
            m = smf.mixedlm(f"{lab} ~ dist_z", dd, groups=dd["patient"].values).fit(reml=True)
            beta = float(m.params["dist_z"])
            ci = m.conf_int().loc["dist_z"].values
            rec.update(beta_per_SD=beta, ci_low=float(ci[0]), ci_high=float(ci[1]),
                       mixedlm_p=float(m.pvalues["dist_z"]))
        except Exception as e:  # noqa: BLE001
            LOG.warning("mixedlm failed for %s (%s)", lab, e)
            rec.update(beta_per_SD=np.nan, ci_low=np.nan, ci_high=np.nan, mixedlm_p=np.nan)
        slopes = [np.polyfit(g["dist_z"], g[lab], 1)[0]
                  for _, g in dd.groupby("patient") if g["dist_z"].nunique() > 2]
        rec["patients_positive_slope"] = int(np.sum(np.array(slopes) > 0))
        rec["patients_tested"] = len(slopes)
        rec["interpretation"] = ("beta>=0 / positive within-patient slopes => NOT enriched toward "
                                 "lesions => inconsistent with paracrine induction")
        a_rows.append(rec)
    a_df = pd.DataFrame(a_rows)
    a_df.to_csv(os.path.join(out, "A_distance_stats_by_geneset.tsv"), sep="\t", index=False)
    LOG.info("(A) distance stats:\n%s", a_df.to_string(index=False))

    # figure: distance decile curves, one line per gene set (z-scored score for comparability)
    fig, ax = plt.subplots(figsize=(7, 5))
    for lab in labels:
        dd = dist_df.dropna(subset=["dist_z", lab]).copy()
        if len(dd) < 50:
            continue
        s = (dd[lab] - dd[lab].mean()) / (dd[lab].std() + 1e-9)
        dd = dd.assign(score_z=s)
        dd["bin"] = pd.qcut(dd["dist_z"].rank(method="first"), 8, labels=False)
        curve = dd.groupby("bin").agg(dist=("dist_z", "mean"), score=("score_z", "mean")).reset_index()
        ax.plot(curve["dist"], curve["score"], marker="o", label=lab)
    ax.axhline(0, color="grey", lw=0.6, ls="--")
    ax.set_xlabel("distance to nearest lesion (within-sample z)")
    ax.set_ylabel("CM-program score (z-scored per set)")
    ax.set_title("(A) distance-to-lesion: flat/positive = cell-autonomous, not paracrine")
    ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "A_distance_to_lesion_by_geneset.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)

    # ---- (B) zone gradient, per gene set ----
    b_rows = []
    for lab in labels:
        sub = spot_df.dropna(subset=[lab])
        groups = {z: sub.loc[sub.zone == z, lab].dropna().values for z in ZONE_ORDER}
        # spot-level (published readout: huge n, kruskal/MWU)
        try:
            _, p_kw = kruskal(*[g for g in groups.values() if len(g)])
        except Exception:  # noqa: BLE001
            p_kw = np.nan

        def _mwu(x, y):
            try:
                return float(mannwhitneyu(x, y, alternative="two-sided")[1])
            except Exception:  # noqa: BLE001
                return np.nan

        # patient-level paired Wilcoxon (matches the reference p~0.001): median per (patient, zone)
        sub2 = sub.copy()
        sub2["patient"] = sub2["sample"].astype(str).str.replace(r"-\d+$", "", regex=True)
        pmed = sub2.groupby(["patient", "zone"])[lab].median().unstack("zone")

        def _paired(zone_a, zone_b):
            if zone_a not in pmed.columns or zone_b not in pmed.columns:
                return np.nan, 0
            pair = pmed[[zone_a, zone_b]].dropna()
            if len(pair) < 3:
                return np.nan, len(pair)
            try:
                return float(wilcoxon(pair[zone_a], pair[zone_b])[1]), len(pair)
            except Exception:  # noqa: BLE001
                return np.nan, len(pair)

        p_pf, n_pf = _paired("preserved", "fibrotic")
        p_pg, n_pg = _paired("preserved", "granulomatous")
        b_rows.append({
            "geneset": lab,
            "mean_preserved": float(np.mean(groups["preserved"])) if len(groups["preserved"]) else np.nan,
            "mean_granulomatous": float(np.mean(groups["granulomatous"])) if len(groups["granulomatous"]) else np.nan,
            "mean_fibrotic": float(np.mean(groups["fibrotic"])) if len(groups["fibrotic"]) else np.nan,
            "median_preserved": float(np.median(groups["preserved"])) if len(groups["preserved"]) else np.nan,
            "median_fibrotic": float(np.median(groups["fibrotic"])) if len(groups["fibrotic"]) else np.nan,
            "median_granulomatous": float(np.median(groups["granulomatous"])) if len(groups["granulomatous"]) else np.nan,
            "kruskal_p_spotlevel": p_kw,
            "mwu_p_preserved_vs_fibrotic_spotlevel": _mwu(groups["preserved"], groups["fibrotic"]),
            "mwu_p_preserved_vs_granuloma_spotlevel": _mwu(groups["preserved"], groups["granulomatous"]),
            "wilcoxon_p_preserved_vs_fibrotic_patientlevel": p_pf,
            "n_patients_pf": n_pf,
            "wilcoxon_p_preserved_vs_granuloma_patientlevel": p_pg,
            "n_patients_pg": n_pg,
        })
    b_df = pd.DataFrame(b_rows)
    b_df.to_csv(os.path.join(out, "B_zone_stats_by_geneset.tsv"), sep="\t", index=False)
    LOG.info("(B) zone stats:\n%s", b_df.to_string(index=False))

    # figure: per-zone mean (z-scored per set) heatmap-ish bar grid
    fig, axes = plt.subplots(1, len(labels), figsize=(3 * len(labels), 4), sharey=True)
    if len(labels) == 1:
        axes = [axes]
    zones_plot = [z for z in ["preserved", "granulomatous", "fibrotic"]]
    for ax, lab in zip(axes, labels):
        sub = spot_df.dropna(subset=[lab])
        s = (sub[lab] - sub[lab].mean()) / (sub[lab].std() + 1e-9)
        sub = sub.assign(score_z=s)
        means = [sub.loc[sub.zone == z, "score_z"].mean() for z in zones_plot]
        ax.bar(range(len(zones_plot)), means, color=["#2c7bb6", "#fdae61", "#d7191c"])
        ax.axhline(0, color="grey", lw=0.6)
        ax.set_xticks(range(len(zones_plot)))
        ax.set_xticklabels(zones_plot, rotation=30, fontsize=7)
        ax.set_title(lab, fontsize=9)
    axes[0].set_ylabel("mean score (z per set)")
    fig.suptitle("(B) zone gradient: preserved-enriched = stage-independent CM program")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "B_zone_gradient_by_geneset.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)

    return a_df, b_df


# ---------------------------------------------------------------------------
# C : snRNA zone flatness on CM pseudobulk (no raw spatial data needed)
# ---------------------------------------------------------------------------
def run_snrna(pb_counts, pb_meta, zones_file, genesets, out):
    import scanpy as sc
    import anndata as ad
    from scipy.stats import kruskal

    counts = pd.read_csv(pb_counts, sep="\t", index_col=0)          # genes x samples
    meta = pd.read_csv(pb_meta, sep="\t", index_col=0)
    zones = pd.read_csv(zones_file, sep="\t")

    counts.index = counts.index.astype(str).str.upper()
    # CS-only samples that have a zone label
    cs_samples = meta.index[meta["disease"].astype(str).str.upper() == "CS"]
    zmap = dict(zip(zones["sample_id"].astype(str), zones["zone"].astype(str)))
    keep = [s for s in counts.columns if s in set(cs_samples) and zmap.get(s, "ambiguous") != "ambiguous"]
    LOG.info("snRNA: %d CS samples with a non-ambiguous zone label", len(keep))

    X = counts[keep].T.to_numpy(dtype=float)                        # samples x genes
    a = ad.AnnData(X=X, obs=pd.DataFrame(index=keep), var=pd.DataFrame(index=counts.index))
    a.obs["zone"] = [zmap[s] for s in keep]
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)

    zorder = ["preserved", "granulomatous", "fibrotic"]
    rows = []
    for lab, genes in genesets:
        present = score_set(a, genes, sc, f"score_{lab}")
        vals = {z: a.obs.loc[a.obs.zone == z, f"score_{lab}"].dropna().values for z in zorder}
        try:
            _, p_kw = kruskal(*[v for v in vals.values() if len(v)])
        except Exception:  # noqa: BLE001
            p_kw = np.nan
        rec = {"geneset": lab, "n_genes_present": len(present),
               "kruskal_p_across_zones": p_kw}
        for z in zorder:
            rec[f"mean_{z}"] = float(np.mean(vals[z])) if len(vals[z]) else np.nan
            rec[f"n_{z}"] = int(len(vals[z]))
        rec["interpretation"] = "kruskal p>0.05 => flat across histologic zones => stage-independent"
        rows.append(rec)
    c_df = pd.DataFrame(rows)
    c_df.to_csv(os.path.join(out, "C_snrna_zone_by_geneset.tsv"), sep="\t", index=False)
    LOG.info("(C) snRNA zone flatness:\n%s", c_df.to_string(index=False))
    return c_df


# ---------------------------------------------------------------------------
# D : cross-disease ordering in CM-dominant spots (identical to cm_spatial_crossdisease.py)
# ---------------------------------------------------------------------------
def run_crossdisease(manifest, genesets, out):
    import scanpy as sc
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.stats import mannwhitneyu

    man = pd.read_csv(manifest, sep="\t", header=None, names=["condition", "path"])
    LOG.info("cross-disease manifest: %s", dict(zip(man.condition, man.path)))
    labels = [lab for lab, _ in genesets]

    rows = []       # per CM-dominant spot: condition, sample, geneset, score
    for _, r in man.iterrows():
        cond, path = str(r.condition), str(r.path)
        try:
            objs = load_geo_visium(path) if os.path.isdir(path) else \
                [(os.path.basename(path).replace(".h5ad", ""), __import__("anndata").read_h5ad(path))]
        except Exception as e:  # noqa: BLE001
            LOG.warning("%s: load failed (%s)", cond, e)
            continue
        for name, a in objs:
            a = ensure_symbols(a)
            a.var_names = a.var_names.astype(str).str.upper()
            a.var_names_make_unique()
            if float(np.asarray(a.X.max())) > 50:
                sc.pp.normalize_total(a, target_sum=1e4)
                sc.pp.log1p(a)
            lin_scores = {}
            for lin, gs in LINEAGES.items():
                gg = [g for g in gs if g in a.var_names]
                if gg:
                    sc.tl.score_genes(a, gg, score_name=f"_lin_{lin}", use_raw=False, random_state=0)
                    lin_scores[lin] = a.obs[f"_lin_{lin}"].values
            if "Cardiomyocyte" not in lin_scores:
                LOG.warning("%s/%s: no CM markers -> skip", cond, name)
                continue
            M = pd.DataFrame(lin_scores)
            cm_dom = (M.idxmax(axis=1).values == "Cardiomyocyte")
            n_cm = int(cm_dom.sum())
            LOG.info("%s/%s: %d/%d CM-dominant spots", cond, name, n_cm, a.n_obs)
            if n_cm == 0:
                continue
            for lab, genes in genesets:
                score_set(a, genes, sc, f"score_{lab}")
                for v in a.obs[f"score_{lab}"].values[cm_dom]:
                    rows.append({"condition": cond, "sample": name, "geneset": lab, "score": float(v)})

    if not rows:
        LOG.warning("no CM-dominant spots scored for cross-disease; skipping D")
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(out, "D_crossdisease_perspot.tsv.gz"), sep="\t", index=False)

    order = [c for c in dict.fromkeys(man.condition.astype(str)) if c in df.condition.unique()]
    # summary: median score per (geneset, condition) + CS-vs-each MWU
    stats = []
    for lab in labels:
        sub = df[df.geneset == lab]
        med = {c: float(np.median(sub.loc[sub.condition == c, "score"])) for c in order
               if len(sub.loc[sub.condition == c])}
        rec = {"geneset": lab}
        for c in order:
            rec[f"median_{c}"] = med.get(c, np.nan)
        if "CS" in order:
            cs = sub.loc[sub.condition == "CS", "score"].values
            for c in order:
                if c == "CS":
                    continue
                other = sub.loc[sub.condition == c, "score"].values
                try:
                    _, p = mannwhitneyu(cs, other, alternative="two-sided")
                except Exception:  # noqa: BLE001
                    p = np.nan
                rec[f"CSvs{c}_dir"] = "CS_higher" if med.get("CS", np.nan) > med.get(c, np.nan) else "CS_lower"
                rec[f"CSvs{c}_p"] = float(p) if p == p else np.nan
        # ordering string (conditions ranked by median, high->low)
        ranked = sorted([c for c in order if c in med], key=lambda c: med[c], reverse=True)
        rec["ordering_high_to_low"] = " > ".join(ranked)
        stats.append(rec)
    d_df = pd.DataFrame(stats)
    d_df.to_csv(os.path.join(out, "D_crossdisease_stats_by_geneset.tsv"), sep="\t", index=False)
    LOG.info("(D) cross-disease:\n%s", d_df.to_string(index=False))

    # figure: median score by condition, one grouped bar cluster per gene set (z per set)
    fig, ax = plt.subplots(figsize=(1.4 * len(order) + 3, 5))
    width = 0.8 / max(1, len(labels))
    xbase = np.arange(len(order))
    for k, lab in enumerate(labels):
        sub = df[df.geneset == lab]
        s = (sub["score"] - sub["score"].mean()) / (sub["score"].std() + 1e-9)
        sub = sub.assign(score_z=s)
        meds = [sub.loc[sub.condition == c, "score_z"].median() if len(sub.loc[sub.condition == c]) else np.nan
                for c in order]
        ax.bar(xbase + k * width, meds, width, label=lab)
    ax.set_xticks(xbase + width * (len(labels) - 1) / 2)
    ax.set_xticklabels(order, rotation=20)
    ax.axhline(0, color="grey", lw=0.6)
    ax.set_ylabel("median CM-program score (z per set)")
    ax.set_title("(D) cross-disease ordering in CM-dominant spots")
    ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "D_crossdisease_by_geneset.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    return d_df


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genesets_dir", default=None, help="dir with geneset_*.txt (default: ../metadata)")
    ap.add_argument("--genesets", nargs="*", default=None, help="explicit geneset file paths (overrides dir)")
    ap.add_argument("--foong_dir", default=None, help="Foong Visium dir (for A/B)")
    ap.add_argument("--pb_counts", default=None, help="CM pseudobulk counts tsv (for C)")
    ap.add_argument("--pb_meta", default=None, help="CM pseudobulk meta tsv (for C)")
    ap.add_argument("--zones", default=None, help="neyazi zones-by-section tsv (for C)")
    ap.add_argument("--manifest", default=None, help="cross-disease manifest tsv (for D)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--analyses", default="ABCD", help="subset of ABCD to run")
    args = ap.parse_args()
    ensure_dir(args.out)

    if args.genesets:
        paths = args.genesets
    else:
        gd = args.genesets_dir or os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "metadata")
        paths = sorted(glob.glob(os.path.join(gd, "geneset_*.txt")))
    genesets = [read_geneset(p) for p in paths]
    if not genesets:
        sys.exit("no gene-set files found")
    LOG.info("gene sets: %s", {lab: len(g) for lab, g in genesets})

    if "A" in args.analyses or "B" in args.analyses:
        if not args.foong_dir:
            sys.exit("--foong_dir required for analyses A/B")
        run_within_cs(args.foong_dir, genesets, args.out)
    if "C" in args.analyses:
        if not (args.pb_counts and args.pb_meta and args.zones):
            sys.exit("--pb_counts/--pb_meta/--zones required for analysis C")
        run_snrna(args.pb_counts, args.pb_meta, args.zones, genesets, args.out)
    if "D" in args.analyses:
        if not args.manifest:
            sys.exit("--manifest required for analysis D")
        run_crossdisease(args.manifest, genesets, args.out)

    LOG.info("done -> %s", args.out)


if __name__ == "__main__":
    main()
