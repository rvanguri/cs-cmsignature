#!/usr/bin/env python3
"""CM-content-normalized per-gene comparison of the CS-CM signature (CS vs control Visium).

Guards against the density critique: TNNI3K/MLIP are cardiomyocyte-enriched, so higher raw expression
in CS cardiomyocyte-dominant spots could just mean those spots are more CM-pure. Here we ask whether
each signature gene is elevated in CS BEYOND what the spot's cardiomyocyte content predicts.

For each dataset we take CM-dominant spots, compute per-spot CM_content = mean expression of canonical
CM structural markers (TNNT2/MYH7/TTN/ACTN2/MYL2/MYBPC3/TNNI3 — disjoint from the signature), then per
signature gene report:
  - raw median (CS vs control)
  - ratio median: expr / CM_content (CS vs control) + Mann-Whitney
  - CM-content-decile-matched effect: within pooled deciles of CM_content, mean(CS) - mean(control),
    averaged over deciles, and how many of 10 deciles favour CS (a density-independent readout)
Plus expr-vs-CM_content scatter with per-group linear fits for the four CM-intrinsic genes.
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import log, ensure_dir  # noqa: E402
from foong_spatial_figures import load_geo_visium  # noqa: E402
from cm_spatial_crossdisease import ensure_symbols, LINEAGES, PERGENE  # noqa: E402
LOG = log("cm_content")

CM_CONTENT = ["TNNT2", "MYH7", "TTN", "ACTN2", "MYL2", "MYBPC3", "TNNI3"]   # disjoint from signature
INTRINSIC = ["GJB7", "TNNI3K", "MLIP", "PANK1"]                             # the four to scatter


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--eps", type=float, default=0.05, help="floor added to CM_content in the ratio")
    args = ap.parse_args()
    ensure_dir(args.out)

    import scanpy as sc
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.stats import mannwhitneyu

    man = pd.read_csv(args.manifest, sep="\t", header=None, names=["condition", "path"])
    recs = []                                   # per CM-dominant spot: condition, CM_content, gene exprs
    for _, r in man.iterrows():
        cond, path = str(r.condition), str(r.path)
        try:
            objs = load_geo_visium(path) if os.path.isdir(path) else \
                [(os.path.basename(path).replace(".h5ad", ""), __import__("anndata").read_h5ad(path))]
        except Exception as e:  # noqa: BLE001
            LOG.warning("%s: load failed (%s)", cond, e); continue
        for name, a in objs:
            a = ensure_symbols(a)
            a.var_names = a.var_names.astype(str).str.upper(); a.var_names_make_unique()
            if float(np.asarray(a.X.max())) > 50:
                sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
            lin = {}
            for L, gs in LINEAGES.items():
                gg = [g for g in gs if g in a.var_names]
                if gg:
                    sc.tl.score_genes(a, gg, score_name=f"_l_{L}", use_raw=False)
                    lin[L] = a.obs[f"_l_{L}"].values
            if "Cardiomyocyte" not in lin:
                LOG.warning("%s/%s: no CM markers -> skip", cond, name); continue
            cm_dom = pd.DataFrame(lin).idxmax(axis=1).values == "Cardiomyocyte"
            if cm_dom.sum() == 0:
                continue
            content_genes = [g for g in CM_CONTENT if g in a.var_names]
            sig_genes = [g for g in PERGENE if g in a.var_names]
            keys = list(dict.fromkeys(content_genes + sig_genes))
            expr = sc.get.obs_df(a[cm_dom], keys=keys)
            cmc = expr[content_genes].mean(axis=1).values
            for g in sig_genes:
                for e, c in zip(expr[g].values, cmc):
                    recs.append({"condition": cond, "sample": name, "gene": g,
                                 "expr": float(e), "cm_content": float(c)})
            LOG.info("%s/%s: %d CM-dominant spots (content genes %d, sig genes %d)",
                     cond, name, int(cm_dom.sum()), len(content_genes), len(sig_genes))

    if not recs:
        sys.exit("no CM-dominant spots collected")
    df = pd.DataFrame(recs)
    df.to_csv(os.path.join(args.out, "cm_content_perspot.tsv.gz"), sep="\t", index=False)

    # ---- per-gene stats: raw, ratio, decile-matched ----
    rows = []
    for g, sub in df.groupby("gene"):
        cs = sub[sub.condition == "CS"]; ct = sub[sub.condition == "control"]
        if len(cs) == 0 or len(ct) == 0:
            continue
        ratio_cs = cs.expr / (cs.cm_content + args.eps)
        ratio_ct = ct.expr / (ct.cm_content + args.eps)
        try:                                    # MWU on per-spot ratios (valid despite zero-inflation)
            _, p_ratio = mannwhitneyu(ratio_cs, ratio_ct, alternative="two-sided")
        except Exception:  # noqa: BLE001
            p_ratio = np.nan
        # CM-content-decile-matched: bin pooled spots, compare MEANS within bin (means handle the
        # zero-inflation of Visium correctly; per-spot medians are ~0 and uninformative here).
        both = sub.copy()
        both["dec"] = pd.qcut(both.cm_content.rank(method="first"), 10, labels=False)
        gm = both.groupby(["dec", "condition"])["expr"].mean().unstack()
        diffs = (gm.get("CS") - gm.get("control")).dropna() if {"CS", "control"} <= set(gm.columns) else pd.Series(dtype=float)
        rows.append({
            "gene": g,
            # MEAN-based (not median): correct for sparse per-spot counts
            "raw_mean_CS": float(cs.expr.mean()), "raw_mean_ctrl": float(ct.expr.mean()),
            "ratio_mean_CS": float(ratio_cs.mean()), "ratio_mean_ctrl": float(ratio_ct.mean()),
            "ratio_p": p_ratio,
            # authoritative direction = the density-controlled, decile-matched effect
            "decile_matched_mean_diff": float(diffs.mean()) if len(diffs) else np.nan,
            "deciles_favoring_CS": int((diffs > 0).sum()) if len(diffs) else 0,
            "adjusted_direction": ("CS_higher" if len(diffs) and diffs.mean() > 0 else "CS_lower_or_ns"),
        })
    stats = pd.DataFrame(rows).sort_values("gene")
    stats.to_csv(os.path.join(args.out, "cm_content_normalized_stats.tsv"), sep="\t", index=False)
    LOG.info("CM-content-normalized per-gene stats:\n%s", stats.to_string(index=False))

    # ---- scatter: expr vs CM_content, CS vs control, with linear fits (the density-control figure) ----
    for g in [x for x in INTRINSIC if x in df.gene.unique()]:
        sub = df[df.gene == g]
        fig, ax = plt.subplots(figsize=(6, 5))
        for cond, col in [("control", "#888888"), ("CS", "#d62728")]:
            s = sub[sub.condition == cond]
            if len(s) == 0:
                continue
            sp = s.sample(min(3000, len(s)), random_state=0)
            ax.scatter(sp.cm_content, sp.expr, s=6, alpha=0.25, c=col, label=cond)
            if s.cm_content.nunique() > 2:
                b, a0 = np.polyfit(s.cm_content, s.expr, 1)
                xs = np.linspace(s.cm_content.min(), s.cm_content.max(), 50)
                ax.plot(xs, b * xs + a0, c=col, lw=2)
        ax.set_xlabel("CM content (mean of CM structural markers)")
        ax.set_ylabel(f"{g} expression"); ax.set_title(f"{g}: expression vs CM content")
        ax.legend()
        fig.tight_layout(); fig.savefig(os.path.join(args.out, f"scatter_{g}.png"), dpi=180,
                                        bbox_inches="tight"); plt.close(fig)
    LOG.info("done -> %s", args.out)


if __name__ == "__main__":
    main()
