#!/usr/bin/env python3
"""Project CS cells onto the HeartMap SCP3689 reference (Plan Step 7c).

scArches/scANVI label-and-embedding TRANSFER of CS cardiomyocytes/macrophages/fibroblasts onto the
HeartMap atlas — NO raw co-embedding (avoids the decontamination mismatch). Reports where CS sits
relative to DCM/HCM/ICM/ARVC/control, and whether CS fibroblasts populate the COL22A1+ (DCM) vs TNC+
(ICM) activated-fibroblast niches.

Requires the HeartMap atlas (SCP3689) present (download_heartmap_scp3689.py). If absent, exits 0
(non-fatal) so the rest of the pipeline still completes.
"""
from __future__ import annotations
import argparse
import glob
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import log, require, ensure_dir, PANELS  # noqa: E402

LOG = log("run_projection")


def find_atlas(atlas_dir: str) -> str | None:
    hits = glob.glob(os.path.join(atlas_dir, "**", "*.h5ad"), recursive=True)
    return hits[0] if hits else None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--integrated", required=True, help="CS-containing integrated h5ad")
    ap.add_argument("--heartmap_atlas", required=True, help="dir with SCP3689 atlas h5ad")
    ap.add_argument("--out", required=True)
    ap.add_argument("--disease_col", default="disease")
    ap.add_argument("--celltype_col", default="predicted_cell_type")
    args = ap.parse_args()
    require(args.integrated, "integrated h5ad")
    ensure_dir(args.out)

    atlas = find_atlas(args.heartmap_atlas)
    if atlas is None:
        LOG.warning("HeartMap atlas not found in %s — skipping projection (non-fatal). "
                    "Run download_heartmap_scp3689.py first.", args.heartmap_atlas)
        sys.exit(0)

    import scanpy as sc
    import scvi  # scArches path uses scvi-tools reference models

    ref = sc.read_h5ad(atlas)
    query = sc.read_h5ad(args.integrated)
    query = query[query.obs[args.disease_col].astype(str).str.upper().eq("CS")].copy()

    # shared features
    common = sorted(set(ref.var_names.str.upper()) & set(query.var_names.str.upper()))
    if len(common) < 100:
        LOG.warning("only %d shared genes between CS query and HeartMap atlas -> skipping projection "
                    "(likely a gene-id mismatch: symbols vs Ensembl. Relabel one to match the other).",
                    len(common))
        sys.exit(0)      # optional step: don't fail the pipeline
    ref.var_names = ref.var_names.str.upper(); query.var_names = query.var_names.str.upper()
    ref, query = ref[:, common].copy(), query[:, common].copy()
    LOG.info("projection on %d shared genes; %d CS query cells", len(common), query.n_obs)

    # TODO(you): if HeartMap ships a trained scANVI/scVI model, load it and use
    # scvi.model.SCANVI.load_query_data(query, ref_model) -> train(weight_decay=0) -> predict().
    # Otherwise train a reference SCANVI on `ref` here, then map the query:
    scvi.model.SCANVI.setup_anndata(ref, batch_key="study" if "study" in ref.obs else None,
                                    labels_key=args.celltype_col if args.celltype_col in ref.obs
                                    else "cell_type", unlabeled_category="Unknown")
    ref_model = scvi.model.SCANVI(ref)
    ref_model.train(max_epochs=200, early_stopping=True)
    q_model = scvi.model.SCANVI.load_query_data(query, ref_model)
    q_model.train(max_epochs=100, weight_decay=0.0)
    query.obs["heartmap_predicted"] = q_model.predict()
    query.obsm["X_scANVI_ref"] = q_model.get_latent_representation(query)

    # composition table: CS cell-state vs each mimic disease label in the reference space
    comp = query.obs["heartmap_predicted"].value_counts(normalize=True).rename("cs_fraction")
    comp.to_csv(os.path.join(args.out, "cs_projection_composition.tsv"), sep="\t")

    # fibroblast niche occupancy (COL22A1 vs TNC)
    niche_rows = []
    fib = query.obs[args.celltype_col].astype(str).str.contains("ibroblast", na=False)
    for niche, genes in PANELS["fib_niche"].items():
        genes = [g for g in genes if g in query.var_names]
        if genes and fib.sum() > 0:
            X = query[fib.values, genes].X
            X = X.toarray() if hasattr(X, "toarray") else X
            niche_rows.append({"niche": niche, "mean_expr_in_CS_fibroblasts": float(X.mean())})
    pd.DataFrame(niche_rows).to_csv(os.path.join(args.out, "cs_fibroblast_niches.tsv"),
                                    sep="\t", index=False)
    query.write_h5ad(os.path.join(args.out, "cs_projected.h5ad"))
    LOG.info("projection complete -> %s", args.out)


if __name__ == "__main__":
    main()
