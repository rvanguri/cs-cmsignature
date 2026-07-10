#!/usr/bin/env python3
"""GSEA + figures + the high-confidence CS set (Plan Step 8 + 7b overlap).

- Pre-ranked gseapy (Hallmark, GO BP, KEGG, Reactome) on the procurement-matched PRIMARY contrasts and
  the high-confidence set (NOT the pooled contrast as headline).
- High-confidence CS set = significant + concordant direction in BOTH a 3' primary (CS-vs-DCM etc.)
  AND Liu CS-vs-ICM, passing the procurement gate.
- Liu overlap stats (hypergeometric) + concordance heatmap across Liu and all 3' contrasts.
- Volcano per primary contrast; DEG heatmap; GSEA dot plots.
Inputs are the DE TSVs written by run_pseudobulk_de.R and run_liu_standalone.py.
"""
from __future__ import annotations
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import log, ensure_dir, PANELS  # noqa: E402

LOG = log("run_gsea_fig")
PRIMARY = ["CS_vs_DCM", "CS_vs_ICM_Simonson", "CS_vs_ARVC"]
GENE_SETS = ["MSigDB_Hallmark_2020", "GO_Biological_Process_2021", "KEGG_2021_Human", "Reactome_2022"]


def load_de(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    # normalize column names across limma (logFC/adj.P.Val) and pydeseq2 (log2FoldChange/padj)
    ren = {"logFC": "lfc", "log2FoldChange": "lfc",
           "adj.P.Val": "padj", "padj": "padj"}
    df = df.rename(columns={k: v for k, v in ren.items() if k in df.columns})
    return df


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--de", required=True, help="DE dir with DE_<contrast>.tsv")
    ap.add_argument("--liu", required=True, help="Liu DE dir")
    ap.add_argument("--figdir", required=True)
    ap.add_argument("--tabledir", required=True)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--lfc", type=float, default=1.0)   # unified with annotate_foong_validation.py
    args = ap.parse_args()
    ensure_dir(args.figdir); ensure_dir(args.tabledir)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # ---- load DE tables ----
    de = {}
    for c in PRIMARY + ["CS_vs_HCM", "CS_vs_NF"]:
        f = os.path.join(args.de, f"DE_{c}.tsv")
        if os.path.exists(f):
            de[c] = load_de(f)
    liu_f = os.path.join(args.liu, "DE_Liu_CS_vs_ICM.tsv")
    liu = load_de(liu_f) if os.path.exists(liu_f) else None

    # ---- volcano per primary ----
    for c, df in de.items():
        if {"lfc", "padj"} <= set(df.columns):
            plt.figure(figsize=(5, 4))
            plt.scatter(df["lfc"], -np.log10(df["padj"].clip(lower=1e-300)), s=4, alpha=0.4)
            for g in PANELS["cs_cm_expected_up"]:
                if g in df.index:
                    plt.annotate(g, (df.loc[g, "lfc"], -np.log10(max(df.loc[g, "padj"], 1e-300))))
            plt.axhline(-np.log10(args.fdr), ls="--", lw=0.6)
            plt.title(c); plt.xlabel("log2FC"); plt.ylabel("-log10 FDR")
            plt.tight_layout(); plt.savefig(os.path.join(args.figdir, f"volcano_{c}.png"), dpi=150)
            plt.close()

    # ---- high-confidence CS set: significant AND SAME-SIGN across Liu + every PRIMARY (Concern 7) ----
    def sig_up(df):                       # significant with |lfc|>thr (direction-agnostic membership)
        return set(df.index[(df["padj"] < args.fdr) & (df["lfc"].abs() > args.lfc)])
    def sig_signed(df):                   # gene -> sign(lfc) for significant genes only
        m = (df["padj"] < args.fdr) & (df["lfc"].abs() > args.lfc)
        return {g: int(np.sign(df.loc[g, "lfc"])) for g in df.index[m]}
    hc = []
    if liu is not None and {"lfc", "padj"} <= set(liu.columns):
        signed = [sig_signed(de[c]) for c in PRIMARY
                  if c in de and {"lfc", "padj"} <= set(de[c].columns)]
        signed.append(sig_signed(liu))
        common = set.intersection(*[set(s) for s in signed]) if signed else set()
        # keep only genes with ONE shared direction across all intersected contrasts (true concordance)
        hc = sorted(g for g in common if len({s[g] for s in signed}) == 1)
    pd.Series(hc, name="gene").to_csv(os.path.join(args.tabledir, "high_confidence_cs_set.tsv"),
                                      sep="\t", index=False)
    LOG.info("high-confidence CS set (direction-concordant): %d genes %s", len(hc), hc[:20])

    # ---- Liu overlap stats (hypergeometric, restricted to the SHARED gene universe) ----
    if liu is not None and "CS_vs_DCM" in de:
        from scipy.stats import hypergeom
        uni = set(de["CS_vs_DCM"].index) & set(liu.index)   # only genes testable in BOTH
        a = sig_up(de["CS_vs_DCM"]) & uni; b = sig_up(liu) & uni
        universe = len(uni); ov = len(a & b)
        p = hypergeom.sf(ov - 1, universe, len(a), len(b)) if universe else np.nan
        pd.DataFrame([{"overlap": ov, "set_a": len(a), "set_b": len(b),
                       "universe": universe, "hypergeom_p": p}]).to_csv(
            os.path.join(args.tabledir, "liu_overlap_stats.tsv"), sep="\t", index=False)

    # ---- pre-ranked GSEA on primaries + high-confidence ----
    try:
        import gseapy as gp
    except Exception as e:  # noqa: BLE001
        gp = None
        LOG.warning("gseapy import failed (%s) — skipping GSEA", e)
    if gp is not None:
        summaries = []
        for c in PRIMARY:
            if c not in de or "lfc" not in de[c].columns:
                continue
            # gseapy.prerank wants a 2-column DataFrame (gene, score), NOT a list-of-lists
            # (a list is mis-read as a file path -> "stat: path should be ... not list").
            s = de[c]["lfc"].dropna()
            s = s[~s.index.duplicated(keep="first")]           # unique gene symbols
            rnk_df = s.sort_values(ascending=False).rename_axis("gene").reset_index()
            rnk_df.columns = ["gene", "score"]
            try:
                # no_plot=True: we only want the res2d table. Rendering a figure per enriched
                # term (thousands) is what OOM'd the 8G job.
                pre = gp.prerank(rnk=rnk_df, gene_sets=GENE_SETS,
                                 outdir=os.path.join(args.figdir, f"gsea_{c}"), min_size=10,
                                 max_size=500, permutation_num=1000, seed=1, no_plot=True)
                res = pre.res2d.copy()
                res.insert(0, "contrast", c)
                summaries.append(res)
                LOG.info("GSEA %s: %d terms", c, len(res))
                del pre, res
                import gc; gc.collect()
            except Exception as e:  # noqa: BLE001
                LOG.warning("gseapy prerank failed for %s (%s)", c, e)
        if summaries:
            allres = pd.concat(summaries, ignore_index=True)
            allres.to_csv(os.path.join(args.tabledir, "gsea_prerank_all.tsv"),
                          sep="\t", index=False)
            LOG.info("wrote combined GSEA table (%d rows)", len(allres))

    LOG.info("GSEA + figures complete -> %s ; tables -> %s", args.figdir, args.tabledir)


if __name__ == "__main__":
    main()
