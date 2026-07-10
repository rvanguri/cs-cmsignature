#!/usr/bin/env python3
"""Memory-light cardiomyocyte pseudobulk from the integrated h5ad (Plan Step 7, aggregation half).

Loading all ~1.1M cells into R (reticulate) OOMs even at 128 GB. Instead we read the h5ad in BACKED
mode, materialize only the CM cells, sum counts per patient in Python, and write a small
(genes x patients) table + patient metadata. run_pseudobulk_de.R then does limma-voom on that tiny
table (no big in-memory object).

Writes:
  <out>/cm_pseudobulk_counts.tsv   genes x patients (integer summed counts)
  <out>/cm_pseudobulk_meta.tsv     one row per patient: disease/study/anatomy/sex + n_cells
"""
from __future__ import annotations
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import log, require, ensure_dir  # noqa: E402

LOG = log("pseudobulk_cm")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--integrated", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--celltype_col", default="predicted_cell_type")
    ap.add_argument("--patient_col", default="individual")
    ap.add_argument("--meta_cols", default="disease,study,anatomy,sex")
    # lineage selection (generalized from CM-only): case-insensitive regex on the cell-type label,
    # and a prefix so CM / T-cell / B-cell runs write to separate files.
    ap.add_argument("--celltype_pattern", default="ardiomyocyte",
                    help="regex matched case-insensitively against the cell-type label")
    ap.add_argument("--label", default="cm",
                    help="output prefix -> <label>_pseudobulk_counts.tsv / _meta.tsv")
    ap.add_argument("--min_cells", type=int, default=10,
                    help="drop a patient's pseudobulk if it has < this many cells of the lineage "
                         "(immune lineages are sparse; tiny pseudobulks are pure noise)")
    ap.add_argument("--subtype_tsv", default=None,
                    help="optional obs_name->refined_label TSV (from subtype_tnk.py); overrides "
                         "the cell-type label for the listed cells before pattern matching")
    args = ap.parse_args()
    require(args.integrated, "integrated h5ad"); ensure_dir(args.out)

    import anndata as ad
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp

    A = ad.read_h5ad(args.integrated, backed="r")
    ct = A.obs[args.celltype_col].astype(str).copy()
    if args.subtype_tsv and os.path.exists(args.subtype_tsv):
        ov = pd.read_csv(args.subtype_tsv, sep="\t", index_col=0)["refined_label"]
        common = ct.index.intersection(ov.index)
        ct.loc[common] = ov.reindex(common).values
        LOG.info("[%s] applied subtype overrides to %d cells from %s",
                 args.label, len(common), os.path.basename(args.subtype_tsv))
    cm = ct.str.contains(args.celltype_pattern, case=False, regex=True, na=False).values
    LOG.info("[%s] cells matching /%s/: %d / %d",
             args.label, args.celltype_pattern, int(cm.sum()), A.n_obs)
    if cm.sum() == 0:
        sys.exit("no cells matched /%s/ in %s (check %s / the label vocabulary)"
                 % (args.celltype_pattern, args.integrated, args.celltype_col))

    sub = A[cm].to_memory()                    # load ONLY the selected-lineage cells
    X = sub.X
    X = sp.csr_matrix(X) if not sp.issparse(X) else X.tocsr()

    pat = sub.obs[args.patient_col].astype(str).values
    cats = pd.unique(pat)
    codes = pd.Categorical(pat, categories=cats).codes
    # indicator (patients x cells) @ (cells x genes) = (patients x genes) summed counts
    ind = sp.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))),
                        shape=(len(cats), len(codes)))
    pb = np.rint(np.asarray((ind @ X).todense())).astype(int)   # patients x genes

    ncell = pd.Series(pat).value_counts().reindex(cats).values
    keep = ncell >= args.min_cells
    if not keep.all():
        LOG.info("[%s] dropping %d/%d patients with < %d cells", args.label,
                 int((~keep).sum()), len(cats), args.min_cells)
    pb, cats, ncell = pb[keep], cats[keep], ncell[keep]
    if len(cats) < 4:
        sys.exit("[%s] only %d patients survive the >=%d-cell filter — too few for DE"
                 % (args.label, len(cats), args.min_cells))

    genes = sub.var_names.astype(str)
    counts = pd.DataFrame(pb.T, index=genes, columns=cats)       # genes x patients
    counts.to_csv(os.path.join(args.out, f"{args.label}_pseudobulk_counts.tsv"), sep="\t")

    mcols = [c for c in args.meta_cols.split(",") if c in sub.obs.columns]
    meta = (sub.obs.drop_duplicates(args.patient_col)
                   .set_index(args.patient_col)[mcols].reindex(cats))
    meta["n_cells"] = ncell
    meta.index.name = args.patient_col
    meta.to_csv(os.path.join(args.out, f"{args.label}_pseudobulk_meta.tsv"), sep="\t")
    LOG.info("[%s] pseudobulk: %d patients x %d genes -> %s",
             args.label, len(cats), len(genes), args.out)


if __name__ == "__main__":
    main()
