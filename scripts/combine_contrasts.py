#!/usr/bin/env python3
"""Combined-evidence table across CS-vs-mimic contrasts for ANY compartment (no external validation).

Same multi-contrast aggregation as the CM/Foong script, but Foong-free: for the T-cell and B-cell
compartments (or any de_<lineage> dir). Ranks genes by how consistently they move across our contrasts,
with an expression term so low-expression fold-change spikes don't dominate.

Per gene:
  n_contrasts_sig  # contrasts where sig (padj<fdr & |lfc|>lfc)
  our_dir          +1/-1 consensus direction across the sig contrasts (0 if mixed)
  our_mean_lfc     mean logFC over the sig contrasts
  mean_aveexpr     mean AveExpr (log2-CPM) across contrasts
  best_contrast/best_lfc/best_padj  the most-significant contrast
  lfc_<contrast>   per-contrast logFC (transparency)
  evidence_score = n_contrasts_sig + expr_bonus,  expr_bonus in [-1,1]
"""
from __future__ import annotations
import argparse, glob, os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import log  # noqa: E402
LOG = log("combine_de")


def load_de(path):
    d = pd.read_csv(path, sep="\t", index_col=0)
    lfc = "logFC" if "logFC" in d.columns else "lfc"
    padj = "adj.P.Val" if "adj.P.Val" in d.columns else "padj"
    cols = {lfc: "lfc", padj: "padj"}
    if "AveExpr" in d.columns:
        cols["AveExpr"] = "aveexpr"
    out = d[list(cols)].rename(columns=cols)
    if "aveexpr" not in out.columns:
        out["aveexpr"] = np.nan
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--de_dir", required=True, help="dir with DE_CS_vs_*.tsv (e.g. de/de_tcell)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--lfc", type=float, default=1.0)
    ap.add_argument("--expr_floor", type=float, default=2.0)
    ap.add_argument("--expr_scale", type=float, default=4.0)
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.de_dir, "DE_CS_vs_*.tsv")))
    files = [f for f in files if "foong" not in os.path.basename(f).lower()]
    if not files:
        sys.exit("no DE_CS_vs_*.tsv in %s" % args.de_dir)
    contrasts = {os.path.basename(f).replace("DE_CS_vs_", "").replace(".tsv", ""): load_de(f) for f in files}
    LOG.info("%s contrasts: %s", args.de_dir, ", ".join(contrasts))

    genes = sorted(set().union(*[d.index for d in contrasts.values()]))
    lfcs = pd.DataFrame({c: d["lfc"].reindex(genes) for c, d in contrasts.items()}, index=genes)
    padjs = pd.DataFrame({c: d["padj"].reindex(genes) for c, d in contrasts.items()}, index=genes)
    mean_aveexpr = pd.DataFrame({c: d["aveexpr"].reindex(genes)
                                 for c, d in contrasts.items()}, index=genes).mean(axis=1)

    sig = (padjs < args.fdr) & (lfcs.abs() > args.lfc)
    n_sig = sig.sum(axis=1)
    sign_in_sig = np.sign(lfcs).where(sig)
    up = (sign_in_sig > 0).sum(axis=1); dn = (sign_in_sig < 0).sum(axis=1)
    our_dir = np.where(up > dn, 1, np.where(dn > up, -1, 0))
    mean_lfc = lfcs.where(sig).mean(axis=1)
    best_c = padjs.idxmin(axis=1); best_padj = padjs.min(axis=1)
    best_lfc = pd.Series([lfcs.loc[g, best_c[g]] for g in genes], index=genes)
    expr_bonus = np.clip((mean_aveexpr - args.expr_floor) / args.expr_scale, -1.0, 1.0).fillna(0)
    ev = n_sig + expr_bonus

    out = pd.DataFrame({
        "n_contrasts_sig": n_sig, "our_dir": our_dir, "our_mean_lfc": mean_lfc.round(3),
        "mean_aveexpr": mean_aveexpr.round(2), "well_expressed": mean_aveexpr >= args.expr_floor,
        "best_contrast": best_c, "best_lfc": best_lfc.round(3), "best_padj": best_padj,
        "evidence_score": ev.round(3),
    }, index=genes)
    for c in contrasts:
        out[f"lfc_{c}"] = lfcs[c].round(3)
    out.index.name = "gene"
    out = out.sort_values(["evidence_score", "best_padj"], ascending=[False, True])
    out.to_csv(args.out, sep="\t")

    LOG.info("genes: %d | contrasts: %d", len(genes), len(contrasts))
    hi = out[(out.n_contrasts_sig >= 2) & (out.our_dir > 0) & out.well_expressed]
    LOG.info("robust up-in-CS (>=2 contrasts, well-expressed): %d", len(hi))
    LOG.info("top 12 by evidence:\n%s",
             out.head(12)[["n_contrasts_sig", "our_dir", "our_mean_lfc", "mean_aveexpr",
                           "evidence_score"]].to_string())
    LOG.info("wrote %s", args.out)


if __name__ == "__main__":
    main()
