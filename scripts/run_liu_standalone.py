#!/usr/bin/env python3
"""Liu CS-vs-ICM 5' co-primary, recomputed FRESH (Plan Step 7b).

The only within-study, within-chemistry, explant-matched CS contrast -> co-primary, NOT co-embedded.
Runs the same QC + (uniform) decontam state + CM pseudobulk + pydeseq2/apeglm shrinkage on Liu's
5' v1.1 data alone, then reports overlap with the 3' primaries (done in run_gsea_figures.py).

This script handles the Liu-only pseudobulk + DE; it assumes Liu has already been decontaminated and
QC'd through the same code paths (run_decontx.R / run_qc.py with dataset='liu'). Pass the QC'd Liu
h5ad via --raw (a dir containing liu_qc.h5ad) or a direct .h5ad.
"""
from __future__ import annotations
import argparse
import os
import sys

# Cap BLAS threads per worker BEFORE numpy import. pydeseq2 parallelizes gene-wise dispersion/LFC
# fitting across workers; nested BLAS threads x workers oversubscribe CPUs on fat nodes and can hang.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import log, require, ensure_dir  # noqa: E402

LOG = log("run_liu")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True, help="dir with liu_qc.h5ad, or a direct .h5ad path")
    ap.add_argument("--out", required=True)
    ap.add_argument("--disease_col", default="disease")
    ap.add_argument("--patient_col", default="individual")
    ap.add_argument("--celltype_col", default="predicted_cell_type")
    ap.add_argument("--n_cpus", type=int, default=None,
                    help="pydeseq2 workers. Default: the SLURM/cgroup allocation. pydeseq2's own "
                         "default (os.cpu_count()) over-spawns on fat nodes and hangs at "
                         "'Fitting dispersions...'; this caps it to what you were actually granted.")
    args = ap.parse_args()
    ensure_dir(args.out)

    import scanpy as sc
    h5 = args.raw if args.raw.endswith(".h5ad") else os.path.join(args.raw, "liu_qc.h5ad")
    require(h5, "Liu QC h5ad")
    adata = sc.read_h5ad(h5)

    # Liu is the standalone 5' arm -> never went through scANVI, so it has no cell-type labels.
    # Annotate by marker-gene scoring (argmax over the lineage panels) to get predicted_cell_type.
    if args.celltype_col not in adata.obs.columns:
        from _common import CELLTYPE_MARKERS
        LOG.info("no '%s' in Liu obs -> marker-based cell typing", args.celltype_col)
        tmp = adata.copy()
        sc.pp.normalize_total(tmp, target_sum=1e4); sc.pp.log1p(tmp)
        sc_cols = {}
        for ct, genes in CELLTYPE_MARKERS.items():
            g = [x for x in genes if x in tmp.var_names]
            if not g:
                continue
            sc.tl.score_genes(tmp, g, score_name="_s")
            sc_cols[ct] = tmp.obs["_s"].values
        if not sc_cols:
            sys.exit("Liu: no marker genes found in var — are var_names gene symbols? "
                     "(run relabel_genes.py on liu_qc.h5ad first)")
        adata.obs[args.celltype_col] = pd.DataFrame(sc_cols, index=adata.obs_names).idxmax(axis=1).values
        LOG.info("Liu cell types: %s", dict(adata.obs[args.celltype_col].value_counts()))

    # cardiomyocytes only
    cm = adata.obs[args.celltype_col].astype(str).str.contains("ardiomyocyte", na=False)
    adata = adata[cm.values].copy()

    # pseudobulk per patient (sum -> integer)
    pb = (
        pd.DataFrame(adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X),
                     index=adata.obs[args.patient_col].astype(str).values,
                     columns=adata.var_names)
        .groupby(level=0).sum().round().astype(int)
    )
    pat_disease = (adata.obs[[args.patient_col, args.disease_col]]
                   .astype(str).drop_duplicates()
                   .set_index(args.patient_col)[args.disease_col])
    pat_disease = pat_disease.reindex(pb.index)

    # pydeseq2 CS vs ICM
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    meta = pd.DataFrame({"condition": pat_disease.str.upper().values}, index=pb.index)
    keep = meta["condition"].isin(["CS", "ICM"])
    # Bound parallelism to the real allocation. sched_getaffinity(0) returns the cgroup-granted CPUs
    # (SLURM --cpus-per-task); os.cpu_count() (pydeseq2's default) sees ALL cores on the node and
    # over-spawns loky workers -> deadlock at "Fitting dispersions...".
    n_cpus = args.n_cpus or len(os.sched_getaffinity(0))
    LOG.info("pydeseq2 n_cpus=%d", n_cpus)
    try:  # pydeseq2 >= 0.4: parallelism configured via an Inference object
        from pydeseq2.default_inference import DefaultInference
        dds = DeseqDataSet(counts=pb.loc[keep], metadata=meta.loc[keep],
                           design_factors="condition", inference=DefaultInference(n_cpus=n_cpus))
    except (ImportError, TypeError):  # older pydeseq2: n_cpus kwarg on DeseqDataSet
        dds = DeseqDataSet(counts=pb.loc[keep], metadata=meta.loc[keep],
                           design_factors="condition", n_cpus=n_cpus)
    dds.deseq2()
    res = DeseqStats(dds, contrast=["condition", "CS", "ICM"])
    res.summary()
    out = os.path.join(args.out, "DE_Liu_CS_vs_ICM.tsv")
    res.results_df.to_csv(out, sep="\t")
    n_sig = int((res.results_df["padj"] < 0.05).sum())
    LOG.info("Liu CS-vs-ICM: %d genes, %d padj<0.05 -> %s", len(res.results_df), n_sig, out)


if __name__ == "__main__":
    main()
