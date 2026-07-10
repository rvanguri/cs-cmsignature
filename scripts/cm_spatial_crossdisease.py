#!/usr/bin/env python3
"""Score the CS cardiomyocyte signature within CM-dominant Visium spots ACROSS diseases.

The CM-core specificity figure: even in cardiomyocyte-dominant regions, does CS express an
inflammasome/arrhythmogenic program that CM regions of other hearts do not? Takes a manifest
(condition<TAB>path; path = a Foong/GEO Visium dir or an .h5ad), and for each dataset:
  1. load/normalize spots
  2. score canonical lineages; keep spots where CARDIOMYOCYTE is the dominant lineage
  3. score the CS-CM up-signature on those CM-dominant spots
Outputs: box/violin of the CS-CM score in CM-dominant spots by condition (with CS-vs-each stats),
per-condition spatial maps of the score restricted to CM-dominant spots, and a summary table.

CS Visium (Foong GSE314910) is ready to use as one condition now; add control (Kuppe) / HCM rows to
the manifest as you stage them.
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import log, ensure_dir  # noqa: E402
from foong_spatial_figures import load_geo_visium, plot_pair  # reuse the Visium loader + paired plot
LOG = log("cm_spatial_xdz")

LINEAGES = {
    "Cardiomyocyte": ["TNNT2", "MYH7", "TTN", "ACTN2", "RYR2", "MYL2", "TNNI3", "MYBPC3"],
    "Fibroblast":    ["DCN", "PDGFRA", "GSN", "COL1A1", "LUM"],
    "Endothelial":   ["PECAM1", "VWF", "CDH5"],
    "Myeloid":       ["CD163", "LYZ", "C1QA", "C1QB"],
    "Lymphoid":      ["CD3E", "MS4A1", "CD79A", "IL7R"],
}
# CANONICAL / SINGLE SOURCE OF TRUTH for the two signature axes (import from here everywhere).
#  - CM_intrinsic: the four data-derived cardiomyocyte genes from the discovery DE (arrhythmogenic/
#    structural, cardiomyocyte-lineage, not expected to ride immune ambient).
#  - "inflammatory": a CURATED granuloma-inflammation REFERENCE program (NOT data-derived) used only
#    as a contrast/positive-control axis to test whether CM_intrinsic co-localizes with inflammation.
#    It mixes inflammasome (NLRC4, IL1RAP) and lymphoid/TLS (BACH2, IL7) markers of the granuloma
#    compartment; it is NOT a cardiomyocyte-intrinsic program. The 4th gene is IL7 (a prior variant in
#    foong_spatial_figures.py used CASP4 — deprecated; IL7 is the version scored in Panel D).
CS_CM_SIGS = {
    "CM_intrinsic":  ["GJB7", "TNNI3K", "MLIP", "PANK1"],
    "inflammatory":  ["NLRC4", "IL1RAP", "BACH2", "IL7"],   # granuloma-inflammation reference (curated)
    "CS_CM_full":    ["GJB7", "TNNI3K", "MLIP", "PANK1", "NLRC4", "IL1RAP", "BACH2", "IL7"],
}
PERGENE = ["GJB7", "TNNI3K", "MLIP", "PANK1", "NLRC4", "IL1RAP", "BACH2", "IL7"]


def ensure_symbols(a):
    """CELLxGENE/Kuppe h5ads index var by Ensembl ID; remap to gene SYMBOLS so marker names match."""
    vn = a.var_names.astype(str)
    if vn.str.startswith("ENSG").mean() > 0.3:
        for col in ("feature_name", "gene_symbols", "gene_symbol", "gene_name", "Symbol", "symbol"):
            if col in a.var.columns:
                a.var["_ensembl"] = a.var_names.astype(str)
                a.var_names = a.var[col].astype(str)
                a.var_names_make_unique()
                LOG.info("  remapped Ensembl -> symbols via var['%s']", col)
                break
    return a


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True, help="TSV: condition<TAB>path (dir or .h5ad)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--signature", default=None, help="optional 1-col gene list; else built-in CS-CM set")
    args = ap.parse_args()
    ensure_dir(args.out)

    import scanpy as sc
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.stats import mannwhitneyu

    sigs = dict(CS_CM_SIGS)
    if args.signature and os.path.exists(args.signature):
        sigs = {"custom": pd.read_csv(args.signature, sep="\t", header=None)[0].astype(str).tolist()}

    man = pd.read_csv(args.manifest, sep="\t", header=None, names=["condition", "path"])
    LOG.info("manifest: %s", dict(zip(man.condition, man.path)))

    rows = []            # per CM-dominant spot: condition, sample, signature, score
    gene_rows = []       # per-gene mean expression in CM-dominant spots, per condition/sample
    for _, r in man.iterrows():
        cond, path = str(r.condition), str(r.path)
        try:
            objs = load_geo_visium(path) if os.path.isdir(path) else \
                [(os.path.basename(path).replace(".h5ad", ""), __import__("anndata").read_h5ad(path))]
        except Exception as e:  # noqa: BLE001
            LOG.warning("%s: load failed (%s)", cond, e); continue
        for name, a in objs:
            a = ensure_symbols(a)                       # Ensembl -> symbol (CELLxGENE/Kuppe)
            a.var_names = a.var_names.astype(str).str.upper(); a.var_names_make_unique()
            if float(np.asarray(a.X.max())) > 50:
                sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
            # lineage scores -> CM-dominant = argmax lineage is Cardiomyocyte
            lin_scores = {}
            for lin, gs in LINEAGES.items():
                gg = [g for g in gs if g in a.var_names]
                if gg:
                    # use_raw=False: CELLxGENE objects keep an Ensembl-indexed .raw that score_genes
                    # would otherwise use, undoing our symbol remap.
                    sc.tl.score_genes(a, gg, score_name=f"_lin_{lin}", use_raw=False)
                    lin_scores[lin] = a.obs[f"_lin_{lin}"].values
            if "Cardiomyocyte" not in lin_scores:
                LOG.warning("%s/%s: no CM markers detected -> skipping", cond, name); continue
            M = pd.DataFrame(lin_scores)
            cm_dom = (M.idxmax(axis=1).values == "Cardiomyocyte")
            n_cm = int(cm_dom.sum())
            present = {s: [g for g in gs if g in a.var_names] for s, gs in sigs.items()}
            LOG.info("%s/%s: %d/%d CM-dominant spots (CM_intrinsic %d, inflammatory %d genes present)",
                     cond, name, n_cm, a.n_obs,
                     len(present.get("CM_intrinsic", [])), len(present.get("inflammatory", [])))
            if n_cm == 0:
                continue
            for signame, sg in present.items():          # score each axis separately
                if not sg:
                    continue
                sc.tl.score_genes(a, sg, score_name="_sig", use_raw=False)
                for v in a.obs["_sig"].values[cm_dom]:
                    rows.append({"condition": cond, "sample": name, "signature": signame, "score": float(v)})
            # per-gene mean in CM-dominant spots (the most interpretable readout)
            import scanpy as _sc
            for g in [x for x in PERGENE if x in a.var_names]:
                col = _sc.get.obs_df(a[cm_dom], keys=[g])[g]
                gene_rows.append({"condition": cond, "sample": name, "gene": g, "mean_expr": float(col.mean())})
            # spatial map: CM-intrinsic score on CM-dominant spots
            if "spatial" in a.obsm and present.get("CM_intrinsic"):
                sub = a[cm_dom].copy()
                sc.tl.score_genes(sub, present["CM_intrinsic"], score_name="CM_intrinsic_score", use_raw=False)
                try:
                    plot_pair(sub, name, "CM_intrinsic_score", f"{cond}:{name} CM-intrinsic (CM spots)",
                              os.path.join(args.out, f"spatial_{cond}_{name}.png"))
                except Exception as e:  # noqa: BLE001
                    LOG.warning("spatial plot %s/%s failed (%s)", cond, name, e)

    if not rows:
        sys.exit("no CM-dominant spots scored — check manifest paths / marker presence")
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(args.out, "cm_dominant_scores.tsv"), sep="\t", index=False)
    genedf = pd.DataFrame(gene_rows)
    genedf.to_csv(os.path.join(args.out, "cm_dominant_pergene.tsv"), sep="\t", index=False)

    order = [c for c in dict.fromkeys(man.condition.astype(str)) if c in df.condition.unique()]
    # ---- one violin panel per signature axis ----
    for signame in df.signature.unique():
        sub = df[df.signature == signame]
        data = [sub.loc[sub.condition == c, "score"].values for c in order]
        fig, ax = plt.subplots(figsize=(1.6 * len(order) + 2, 5))
        ax.violinplot(data, showmeans=True, showextrema=False)
        ax.boxplot(data, widths=0.15, showfliers=False)
        ax.set_xticks(range(1, len(order) + 1)); ax.set_xticklabels(order, rotation=20)
        ax.set_ylabel(f"{signame} score"); ax.set_title(f"{signame} in CM-dominant Visium spots")
        fig.tight_layout(); fig.savefig(os.path.join(args.out, f"score_{signame}.png"), dpi=200,
                                        bbox_inches="tight"); plt.close(fig)

    # ---- stats per axis: CS vs each condition, TWO-SIDED (report direction honestly) ----
    stats = []
    if "CS" in order:
        for signame in df.signature.unique():
            sub = df[df.signature == signame]
            cs = sub.loc[sub.condition == "CS", "score"].values
            for c in order:
                if c == "CS":
                    continue
                other = sub.loc[sub.condition == c, "score"].values
                try:
                    _, p = mannwhitneyu(cs, other, alternative="two-sided")
                except Exception:  # noqa: BLE001
                    p = np.nan
                stats.append({"signature": signame, "vs": c,
                              "median_CS": float(np.median(cs)), "median_other": float(np.median(other)),
                              "direction": "CS_higher" if np.median(cs) > np.median(other) else "CS_lower",
                              "mannwhitney_p": p})
    pd.DataFrame(stats).to_csv(os.path.join(args.out, "cm_score_stats.tsv"), sep="\t", index=False)
    LOG.info("median score by signature x condition:\n%s",
             df.groupby(["signature", "condition"])["score"].median().round(3).to_string())
    if len(genedf):
        LOG.info("per-gene mean in CM-dominant spots (gene x condition):\n%s",
                 genedf.groupby(["gene", "condition"])["mean_expr"].mean().round(3).unstack().to_string())
    LOG.info("done -> %s", args.out)


if __name__ == "__main__":
    main()
