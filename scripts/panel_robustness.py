#!/usr/bin/env python3
"""Robustness of the four-gene CS-CM panel to scoring/threshold choices.

The evidence ranking uses hand-set thresholds (fdr, lfc, expression floor). This checks that the
canonical panel (GJB7, TNNI3K, MLIP, PANK1) is stable across a grid of reasonable settings, and gives
a permutation null for the Foong spatial recovery (so panel membership is not an artifact of tuning).
Uses only local DE tables + the Foong significant-DEG csv.
"""
import argparse, glob, os, sys, itertools
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import log  # noqa: E402
LOG = log("panel_robust")
PANEL = ["GJB7", "TNNI3K", "MLIP", "PANK1"]


def load(de_dir):
    con = {}
    for f in glob.glob(os.path.join(de_dir, "DE_CS_vs_*.tsv")):
        if "foong" in f.lower():
            continue
        c = os.path.basename(f).replace("DE_CS_vs_", "").replace(".tsv", "")
        d = pd.read_csv(f, sep="\t", index_col=0)
        con[c] = d.rename(columns={"logFC": "lfc", "adj.P.Val": "padj", "AveExpr": "ae"})[["lfc", "padj", "ae"]]
    return con


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--de_dir", default="results/de/de_cm")
    ap.add_argument("--foong", default=os.environ.get("FOONG_SIGDEG", "SIGDEG.csv"))
    ap.add_argument("--out", default="results/de/de_cm/panel_robustness.tsv")
    ap.add_argument("--nperm", type=int, default=10000)
    args = ap.parse_args()

    con = load(args.de_dir)
    genes = sorted(set().union(*[d.index for d in con.values()]))
    lfcs = pd.DataFrame({c: d["lfc"].reindex(genes) for c, d in con.items()}, index=genes)
    padjs = pd.DataFrame({c: d["padj"].reindex(genes) for c, d in con.items()}, index=genes)
    ae = pd.DataFrame({c: d["ae"].reindex(genes) for c, d in con.items()}, index=genes).mean(1)
    fo = pd.read_csv(args.foong).set_index("gene")
    fo_up = fo["log2FoldChange"].reindex(genes) > 0
    fo_sig = fo["padj"].reindex(genes).notna()          # Foong file is significant-only
    fo_repl_up = (fo_sig & fo_up).fillna(False)

    # ---- grid stability: fraction of settings where each panel gene is 'robust up-in-CS' ----
    grid = list(itertools.product([0.01, 0.05, 0.10], [0.5, 1.0, 1.5], [1.0, 2.0, 3.0]))
    rows = []
    robust_sets = []
    for fdr, lfc, floor in grid:
        up_sig = ((padjs < fdr) & (lfcs > lfc)).sum(1)          # # contrasts up-significant
        robust = (up_sig >= 2) & (ae >= floor) & fo_repl_up      # >=2 concordant-up + well-expr + Foong-up
        robust_sets.append(set(np.array(genes)[robust.values]))
    n = len(grid)
    for g in PANEL:
        frac = np.mean([g in s for s in robust_sets])
        rows.append({"gene": g, "settings_robust": f"{int(frac*n)}/{n}", "fraction": round(frac, 3),
                     "mean_aveexpr": round(float(ae.get(g, np.nan)), 2),
                     "foong_up_replicated": bool(fo_repl_up.get(g, False))})
    stab = pd.DataFrame(rows)

    # ---- null: is 'all 4 panel genes recovered up in Foong' more than chance? ----
    # This tests directional concordance of the FIXED biological panel in INDEPENDENT spatial data
    # (Foong). It does NOT re-draw from the evidence-score ranking that could select the panel (the
    # panel is fixed a priori), and the base rate is computed over up-candidates EXCLUDING the panel
    # genes, so the reference is not inflated by the four genes' own Foong status.
    up_any = ((padjs < 0.05) & (lfcs > 1.0)).sum(1) >= 2       # up-in-CS candidate universe (our data)
    cand = [g for g in np.array(genes)[up_any.values] if g not in PANEL]   # exclude the panel itself
    base_rate = float(fo_repl_up.reindex(cand).mean())        # P(a non-panel up-candidate is Foong-up)
    obs = int(fo_repl_up.reindex(PANEL).sum())
    rng = np.random.default_rng(0)
    # draw 4 genes at random from the (panel-excluded) up-candidate universe, count Foong-up
    fo_arr = fo_repl_up.reindex(cand).fillna(False).values
    null = np.array([rng.choice(fo_arr, size=len(PANEL), replace=False).sum() for _ in range(args.nperm)])
    p_perm = float((null >= obs).mean())
    LOG.info("grid stability (panel genes robust across %d settings):\n%s", n, stab.to_string(index=False))
    LOG.info("permutation null: base Foong-up rate among up-candidates = %.2f; "
             "%d/4 panel genes Foong-up; P(>=%d by chance) = %.4f", base_rate, obs, obs, p_perm)
    stab.to_csv(args.out, sep="\t", index=False)
    with open(args.out.replace(".tsv", "_null.txt"), "w") as fh:
        fh.write(f"base_foong_up_rate={base_rate:.3f} obs={obs}/4 p_perm={p_perm:.4f} nperm={args.nperm}\n")
    print("wrote", args.out)


if __name__ == "__main__":
    main()
