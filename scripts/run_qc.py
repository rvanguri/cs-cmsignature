#!/usr/bin/env python3
"""QC + doublet detection on decontaminated counts (Plan Step 2c), one dataset per call.

Thresholds verbatim from the Plan:
  min_genes >= 200, max_genes <= 8000, MT% < 5, (CellBender cell-prob > 0.5 if present),
  doublets via scrublet on decontaminated counts, keep genes in >= 3 nuclei/sample,
  biotype = protein_coding + lncRNA (applied at gene-harmonization, Step 4).
Writes a QC'd .h5ad.
"""
from __future__ import annotations
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import log, require, ensure_dir, normalize_disease, CANONICAL_DISEASES  # noqa: E402

LOG = log("run_qc")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--counts", required=True, help="DecontX output .h5ad for this dataset")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min_genes", type=int, default=200)
    ap.add_argument("--max_genes", type=int, default=8000)
    ap.add_argument("--max_mt", type=float, default=5.0)
    ap.add_argument("--sample_meta", default="",
                    help="sample_meta.tsv to merge disease/procurement/region/sex/individual onto obs "
                         "(keyed by sample_id; non-destructive: won't overwrite columns already present)")
    args = ap.parse_args()
    require(args.counts, "decontaminated counts")
    ensure_dir(os.path.dirname(args.out) or ".")

    import pandas as pd
    import scanpy as sc
    adata = sc.read_h5ad(args.counts)
    adata.obs["dataset"] = args.dataset
    adata.obs["study"] = args.dataset          # DE 'study' covariate (scANVI also sets this on concat)
    n0 = adata.n_obs

    # --- merge sample-level metadata onto obs (keyed by sample_id from DecontX) ---
    # Attaches disease/anatomy/sex/individual for the GEO + Liu datasets, which carry no embedded
    # metadata. Non-destructive: a column already present (e.g. Reichart's own disease/region from
    # CELLxGENE) is left untouched. region -> anatomy to match the DE model.
    if args.sample_meta and os.path.exists(args.sample_meta) and "sample_id" in adata.obs:
        meta = (pd.read_csv(args.sample_meta, sep="\t", comment="#")
                  .drop_duplicates("sample_id").set_index("sample_id"))
        rename = {"region": "anatomy"}
        for col in ("disease", "procurement", "region", "sex", "individual"):
            if col not in meta.columns:
                continue
            target = rename.get(col, col)
            if target in adata.obs.columns:        # don't clobber metadata already present
                continue
            adata.obs[target] = adata.obs["sample_id"].map(meta[col]).values
            LOG.info("merged '%s' -> obs['%s']", col, target)
        # individual defaults to sample_id when absent/unmapped
        if "individual" not in adata.obs.columns:
            adata.obs["individual"] = adata.obs["sample_id"].values
        else:
            # cast to object first: a dataset carrying 'individual' as a pandas Categorical (e.g. Reichart
            # from CELLxGENE) raises "Cannot set a Categorical ... without identical categories" when
            # filling NaNs with sample_id values outside its category set. Fill on plain strings instead.
            adata.obs["individual"] = (adata.obs["individual"].astype(object)
                                       .fillna(adata.obs["sample_id"].astype(object)).astype(str))
    elif "sample_id" not in adata.obs:
        LOG.warning("no 'sample_id' in obs: skipping sample_meta merge (Reichart/atlas carry own obs)")

    # --- harmonize disease labels to one vocab (CELLxGENE ontology strings -> DCM/ARVC/NF/...) ---
    if "disease" in adata.obs.columns:
        adata.obs["disease"] = adata.obs["disease"].astype(str).map(normalize_disease)
        bad = sorted(set(adata.obs["disease"]) - CANONICAL_DISEASES)
        if bad:
            LOG.warning("%s: disease values not in canonical vocab %s -> %s "
                        "(add aliases in _common.DISEASE_ALIASES)", args.dataset,
                        sorted(CANONICAL_DISEASES), bad)

    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)

    # per-nucleus filters
    sc.pp.filter_cells(adata, min_genes=args.min_genes)
    adata = adata[adata.obs["n_genes_by_counts"] <= args.max_genes].copy()
    adata = adata[adata.obs["pct_counts_mt"] < args.max_mt].copy()

    # CellBender cell-probability gate if the column was carried through
    for col in ("cell_probability", "cellbender_cell_prob"):
        if col in adata.obs:
            adata = adata[adata.obs[col] > 0.5].copy()
            break

    # gene filter: present in >= 3 nuclei
    sc.pp.filter_genes(adata, min_cells=3)

    # doublets on decontaminated counts. Scrublet simulates doublets by DOUBLING the matrix, so its peak
    # memory scales with the largest object it processes. Run it PER SAMPLE (batch_key="sample_id") so the
    # peak is bounded by the largest single sample, not the whole cohort: a multi-sample dataset like
    # neyazi (41 samples, ~283k nuclei total but only ~7-15k per sample) then gets doublet-filtered exactly
    # as intended without OOM. A pre-curated single mega-object with no per-sample structure (e.g. Reichart,
    # one ~880k-nucleus CELLxGENE object, no sample_id) has no batching to bound the peak, so it is skipped
    # (already doublet-filtered upstream). The guard is therefore on the largest BATCH, not total cells.
    SCRUBLET_MAX_BATCH = 200_000
    batch_key = "sample_id" if ("sample_id" in adata.obs and adata.obs["sample_id"].nunique() > 1) else None
    max_batch = int(adata.obs[batch_key].value_counts().max()) if batch_key else adata.n_obs
    if max_batch > SCRUBLET_MAX_BATCH:
        LOG.warning("%s: largest scrublet batch %d > %d -> skipping scrublet "
                    "(pre-curated mega-object; treat as already doublet-filtered)",
                    args.dataset, max_batch, SCRUBLET_MAX_BATCH)
    else:
        try:
            sc.external.pp.scrublet(adata, batch_key=batch_key)  # per-sample; adds 'predicted_doublet'
            if "predicted_doublet" in adata.obs:
                n_before = adata.n_obs
                adata = adata[~adata.obs["predicted_doublet"].astype(bool)].copy()
                LOG.info("%s: scrublet removed %d doublets (batch_key=%s, %d batches)",
                         args.dataset, n_before - adata.n_obs, batch_key,
                         adata.obs[batch_key].nunique() if batch_key else 1)
        except Exception as e:  # noqa: BLE001
            LOG.warning("scrublet skipped (%s): run scDblFinder in R as fallback", e)

    adata.write_h5ad(args.out)
    LOG.info("%s QC: %d -> %d nuclei -> %s", args.dataset, n0, adata.n_obs, args.out)


if __name__ == "__main__":
    main()
