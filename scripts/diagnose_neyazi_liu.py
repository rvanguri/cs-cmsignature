#!/usr/bin/env python3
"""Why do Neyazi (3') and Liu (5') CS DEGs barely overlap? Diagnostic (not part of the DE chain).

Loads the integrated h5ad (Neyazi = study 'neyazi') and liu_qc.h5ad, and reports, per lineage:
  - donors per disease (the real driver of pseudobulk power)
  - cells per disease
  - sequencing depth (median counts & genes per cell) 3' vs 5'
  - detection rate of a few CS-CM candidate genes
so we can see whether low overlap is power (few donors/cells), depth (5' shallower), or coverage.
"""
from __future__ import annotations
import argparse, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import log  # noqa: E402
LOG = log("diagnose")

CANDIDATES = ["TNNI3K", "GJB7", "MLIP", "PANK1", "SFXN4"]


def summarize(A, name, disease_col, ct_col, patient_col):
    import numpy as np, pandas as pd, scipy.sparse as sp
    print(f"\n########## {name}  (n_obs={A.n_obs}, n_vars={A.n_vars}) ##########")
    obs = A.obs
    for c in (disease_col, ct_col, patient_col):
        print(f"  {c}: present={c in obs.columns}")
    # depth
    X = A.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    counts = np.asarray(X.sum(1)).ravel()
    ngenes = np.asarray((X > 0).sum(1)).ravel()
    print(f"  depth: median UMIs/cell={np.median(counts):.0f}  median genes/cell={np.median(ngenes):.0f}")
    if disease_col in obs:
        print("  donors per disease:")
        for d, sub in obs.groupby(disease_col):
            nd = sub[patient_col].nunique() if patient_col in obs else float("nan")
            print(f"     {str(d):12s} donors={nd}  cells={len(sub)}")
    if ct_col in obs:
        print("  cells per lineage:", dict(obs[ct_col].astype(str).value_counts().head(8)))
    # detection of candidate genes among cardiomyocytes
    vn = A.var_names.astype(str)
    for g in CANDIDATES:
        if g in set(vn):
            j = list(vn).index(g)
            col = np.asarray(X[:, j].todense()).ravel()
            print(f"  {g}: detected in {100*(col>0).mean():.1f}% of cells, mean count {col.mean():.3f}")
        else:
            print(f"  {g}: ABSENT from var_names")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--integrated", required=True)
    ap.add_argument("--liu", required=True)
    ap.add_argument("--disease_col", default="disease")
    ap.add_argument("--celltype_col", default="predicted_cell_type")
    ap.add_argument("--patient_col", default="individual")
    args = ap.parse_args()
    import anndata as ad
    A = ad.read_h5ad(args.integrated, backed="r")
    ney = A[A.obs["study"].astype(str).eq("neyazi")].to_memory() if "study" in A.obs else A.to_memory()
    summarize(ney, "NEYAZI (3', from integrated)", args.disease_col, args.celltype_col, args.patient_col)
    L = ad.read_h5ad(args.liu)
    # Liu patient/celltype cols may differ; try a few
    pcol = next((c for c in [args.patient_col, "individual", "sample_id", "donor", "patient"] if c in L.obs), args.patient_col)
    ccol = next((c for c in [args.celltype_col, "predicted_cell_type", "cell_type"] if c in L.obs), args.celltype_col)
    summarize(L, "LIU (5')", args.disease_col, ccol, pcol)


if __name__ == "__main__":
    main()
