#!/usr/bin/env python3
"""Gene-ID harmonization + feature intersection + scANVI integration of the 3' arm (Plan Steps 3-5).

- Concatenate the five QC'd 3' datasets (Neyazi, Larson, Chin2022, Chin2021, Reichart).
- Harmonize to HGNC symbols, intersect features (~15-18k genes).
- Reichart Cell-Ontology labels as the annotation reference; others = "Unknown".
- scANVI: batch_key=study, labels_key=cell_type_harmonized, train(max_epochs=400, early_stopping).
- Liu (5') is NOT here (standalone, run_liu_standalone.py).
Writes the trained model dir + integrated .h5ad (with predicted_cell_type + X_scANVI).
"""
from __future__ import annotations
import argparse
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import log, require, ensure_dir, THREE_PRIME_DATASETS  # noqa: E402

LOG = log("run_scanvi")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--qc_dir", required=True, help="dir with <dataset>_qc.h5ad files")
    ap.add_argument("--reichart_labels", default="reichart",
                    help="dataset name whose obs carries reference cell-type labels")
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out_model", required=True)
    ap.add_argument("--out_h5ad", required=True)
    ap.add_argument("--max_epochs", type=int, default=0,
                    help="0/auto = scale to dataset size (400 epochs is too many for ~1M cells)")
    args = ap.parse_args()
    ensure_dir(args.workdir)
    ensure_dir(os.path.dirname(args.out_h5ad) or ".")

    import anndata as ad
    import scanpy as sc
    import scvi

    # ---- load + concat the five 3' datasets ----
    adatas = {}
    for ds in THREE_PRIME_DATASETS:
        f = os.path.join(args.qc_dir, f"{ds}_qc.h5ad")
        require(f, f"{ds} QC h5ad")
        a = sc.read_h5ad(f)
        a.var_names = a.var_names.str.upper()           # crude HGNC harmonization
        a.var_names_make_unique()
        a.obs["study"] = ds
        adatas[ds] = a

    # feature intersection across all five
    common = set(adatas[THREE_PRIME_DATASETS[0]].var_names)
    for ds in THREE_PRIME_DATASETS[1:]:
        common &= set(adatas[ds].var_names)
    common = sorted(common)
    LOG.info("feature intersection: %d genes", len(common))
    adata = ad.concat({ds: a[:, common] for ds, a in adatas.items()},
                      join="outer", label="study_concat", index_unique="-")
    del adatas

    # Restricting to the common gene set can zero-out low-count cells (their reads were in genes not
    # shared across datasets). scVI's encoder produces NaN on zero-library cells and training diverges,
    # so drop cells with too few counts/genes in the shared space BEFORE training.
    n_before = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=3)
    sc.pp.filter_cells(adata, min_counts=1)
    if adata.n_obs < n_before:
        LOG.info("dropped %d cells with ~0 counts after gene intersection (%d -> %d)",
                 n_before - adata.n_obs, n_before, adata.n_obs)

    # ---- labels: reference dataset provides them (normalized to the shared vocab), others "Unknown" ----
    from _common import normalize_celltype  # noqa: E402
    lab_col = None
    for cand in ("cell_type", "cell_type_harmonized", "Cell_type"):
        if cand in adata.obs:
            lab_col = cand
            break
    adata.obs["cell_type_harmonized"] = "Unknown"
    if lab_col:
        ref = adata.obs["study"] == args.reichart_labels
        adata.obs.loc[ref, "cell_type_harmonized"] = (
            adata.obs.loc[ref, lab_col].astype(str).map(normalize_celltype))
    n_labeled = int((adata.obs["cell_type_harmonized"] != "Unknown").sum())
    LOG.info("reference-labeled cells: %d (%d cell types)", n_labeled,
             adata.obs.loc[adata.obs["cell_type_harmonized"] != "Unknown",
                           "cell_type_harmonized"].nunique())
    if n_labeled == 0:
        sys.exit("ERROR: no reference-labeled cells — scANVI needs >=1 labeled category. "
                 "Check that Reichart's decontx h5ad carries a 'cell_type' column and study=='%s'."
                 % args.reichart_labels)

    # ---- scANVI ----
    scvi.model.SCANVI.setup_anndata(adata, batch_key="study",
                                    labels_key="cell_type_harmonized",
                                    unlabeled_category="Unknown")
    model = scvi.model.SCANVI(adata)
    # Optional dataloader parallelism (default OFF). The bottleneck is GPU compute, not data loading,
    # and num_workers>0 uses /dev/shm which can be capped on HPC nodes (-> Bus error/hang). Opt in with
    # e.g. `export SCANVI_NUM_WORKERS=7` if you find data loading is starving the GPU.
    _nw = int(os.environ.get("SCANVI_NUM_WORKERS", "0"))
    if _nw > 0:
        scvi.settings.dl_num_workers = _nw
        LOG.info("using %d dataloader workers", _nw)
    # Scale epochs to dataset size: 400 epochs on ~1M cells is ~33h (times out). More cells => fewer
    # epochs needed. --max_epochs<=0 means "auto". Gradient clipping guards against NaN divergence.
    import numpy as np
    if args.max_epochs and args.max_epochs > 0:
        epochs = args.max_epochs
    else:
        epochs = int(np.clip(round(1e5 / max(adata.n_obs, 1) * 400), 40, 400))
    LOG.info("training scANVI for %d epochs on %d nuclei", epochs, adata.n_obs)
    try:
        model.train(max_epochs=epochs, early_stopping=True,
                    gradient_clip_val=1.0, plan_kwargs={"lr": 1e-3})
    except TypeError as e:  # older scvi may not accept gradient_clip_val as a trainer kwarg
        LOG.warning("train kwargs not all accepted (%s) -> retrying without gradient clipping", e)
        model.train(max_epochs=epochs, early_stopping=True, plan_kwargs={"lr": 1e-3})
    adata.obsm["X_scANVI"] = model.get_latent_representation()
    adata.obs["predicted_cell_type"] = model.predict()

    model.save(args.out_model, overwrite=True)
    # anndata refuses to write when the obs index name collides with a column of different values
    # (happens after concat suffixes the index). Clear the index name to avoid the clash.
    if adata.obs.index.name is not None and adata.obs.index.name in adata.obs.columns:
        adata.obs = adata.obs.rename_axis(None)
    adata.write_h5ad(args.out_h5ad)
    LOG.info("scANVI done: %d nuclei -> %s ; model -> %s",
             adata.n_obs, args.out_h5ad, args.out_model)


if __name__ == "__main__":
    main()
