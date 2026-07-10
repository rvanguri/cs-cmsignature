#!/usr/bin/env python3
"""Split the scANVI 'TNK' (lymphocyte) compartment into T/NK vs B cells by marker score.

scANVI predicts a single lumped lymphocyte label ('TNK') because the Reichart reference folded
B cells into a generic lymphocyte class (the `lymphocyte -> TNK` alias). Pseudobulking that mixed
bin conflates composition with expression (e.g. MS4A1 looked 'up in CS' only because the T/NK
fraction shifted). This script re-types just the TNK cells: score a B-cell signature vs a T/NK
signature and assign each cell to 'Bcell' or 'Tcell'. Output is a small TSV (obs_name -> refined
label) that pseudobulk_cm.py applies via --subtype_tsv, so we get a clean Tcell and a clean Bcell
compartment without rewriting the 1.1M-cell h5ad.
"""
from __future__ import annotations
import argparse, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import log, require, ensure_dir  # noqa: E402
LOG = log("subtype_tnk")

# B/plasma vs T/NK signatures (kept broad so scoring is robust to dropout)
B_MARKERS = ["MS4A1", "CD79A", "CD79B", "EBF1", "BANK1", "MZB1", "IGHM", "IGKC", "IGHG1", "TNFRSF13B"]
T_MARKERS = ["CD3E", "CD3D", "CD2", "IL7R", "THEMIS", "TRAC", "TRBC2", "CCL5", "NKG7", "GZMK",
             "GZMB", "KLRD1", "CD8A", "CD247"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--integrated", required=True)
    ap.add_argument("--out", required=True, help="output TSV: obs_name<TAB>refined_label")
    ap.add_argument("--celltype_col", default="predicted_cell_type")
    ap.add_argument("--disease_col", default="disease")
    ap.add_argument("--tnk_pattern", default="TNK|T cell|NK|lymphocyte")
    args = ap.parse_args()
    require(args.integrated, "integrated h5ad"); ensure_dir(os.path.dirname(args.out) or ".")

    import anndata as ad
    import scanpy as sc
    import numpy as np
    import pandas as pd

    A = ad.read_h5ad(args.integrated, backed="r")
    mask = A.obs[args.celltype_col].astype(str).str.contains(
        args.tnk_pattern, case=False, regex=True, na=False).values
    LOG.info("TNK/lymphocyte cells: %d / %d", int(mask.sum()), A.n_obs)
    if mask.sum() == 0:
        sys.exit("no TNK/lymphocyte cells matched /%s/" % args.tnk_pattern)

    sub = A[mask].to_memory()
    # X is raw counts (pseudobulk sums+rounds it) -> normalize+log for scoring
    sc.pp.normalize_total(sub, target_sum=1e4)
    sc.pp.log1p(sub)
    Bp = [g for g in B_MARKERS if g in sub.var_names]
    Tp = [g for g in T_MARKERS if g in sub.var_names]
    LOG.info("scoring with %d B markers, %d T/NK markers", len(Bp), len(Tp))
    sc.tl.score_genes(sub, Bp, score_name="B_score", ctrl_size=50)
    sc.tl.score_genes(sub, Tp, score_name="T_score", ctrl_size=50)

    refined = np.where(sub.obs["B_score"].values > sub.obs["T_score"].values, "Bcell", "Tcell")
    out = pd.DataFrame({"refined_label": refined}, index=sub.obs_names)
    out.index.name = "obs_name"
    out.to_csv(args.out, sep="\t")
    LOG.info("wrote %s (%d cells)", args.out, len(out))

    # composition report (this is itself informative: B fraction of lymphocytes per disease)
    if args.disease_col in sub.obs:
        comp = (pd.crosstab(sub.obs[args.disease_col].astype(str), refined, normalize="index") * 100).round(1)
        LOG.info("B/T composition of the lymphocyte compartment by disease (%%):\n%s", comp.to_string())


if __name__ == "__main__":
    main()
