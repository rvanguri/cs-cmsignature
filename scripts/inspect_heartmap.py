#!/usr/bin/env python3
"""Inspect the 80 GB HeartMap_V1.0.h5ad (SCP3689) WITHOUT loading the count matrices, and optionally
write a slim reference for the projection step (Plan Step 7c).

The full file is ~80 GB (2.4M nuclei, 9 studies, dual CellBender+CellRanger count layers). Loading it
whole needs a high-mem node. This opens it in backed mode (metadata only) so you can:
  - confirm it's the integrated object and not a corrupt/partial download,
  - see which obsm embedding to project CS onto (e.g. X_scANVI),
  - read the obs label columns (cell type / disease / study / anatomy) the gates + DE filters expect.

Then `--slim` writes a few-GB reference carrying obs + var + obsm (+ optionally ONE counts layer) so
run_heartmap_projection.py never has to touch the 80 GB file.

Usage:
  # inspect only (seconds; safe on a login/cpu_short node)
  python inspect_heartmap.py --h5ad $PROJ/heartmap/atlas/HeartMap_V1.0.h5ad

  # write slim reference for projection (run on a high-mem fn_* node if --keep-layer is used)
  python inspect_heartmap.py --h5ad .../HeartMap_V1.0.h5ad --slim $PROJ/heartmap/atlas/heartmap_slim.h5ad
"""
from __future__ import annotations
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import log, require  # noqa: E402

LOG = log("inspect_heartmap")

# obs columns we care about downstream; we print value_counts for any that are present
LABEL_CANDIDATES = [
    "cell_type", "cell_type_harmonized", "CellType", "cell_states", "cell_state",
    "disease", "disease_state", "condition",
    "study", "dataset", "anatomy", "region", "anatomical_region",
    "sex", "individual", "donor_id", "donor",
]


def human(n: int) -> str:
    for unit in ("", "K", "M", "B"):
        if abs(n) < 1000:
            return f"{n:.0f}{unit}"
        n /= 1000.0
    return f"{n:.1f}T"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True, help="path to HeartMap_V1.0.h5ad")
    ap.add_argument("--slim", default="", help="if set, write a slim reference h5ad here")
    ap.add_argument("--keep-layer", default="",
                    help="name of ONE counts layer to retain in the slim ref (default: drop all "
                         "matrices, keep obs+var+obsm only). Requires a high-mem node.")
    ap.add_argument("--head", type=int, default=20, help="max categories to print per obs column")
    args = ap.parse_args()
    require(args.h5ad, "HeartMap h5ad")

    import anndata as ad
    import h5py

    # ---- 1) structural peek via h5py (no AnnData parsing of the big matrices) ----
    LOG.info("file: %s (%.1f GB on disk)", args.h5ad, os.path.getsize(args.h5ad) / 1e9)
    with h5py.File(args.h5ad, "r") as f:
        top = list(f.keys())
        LOG.info("top-level groups: %s", top)
        if "layers" in f:
            LOG.info("layers: %s", list(f["layers"].keys()))
        if "obsm" in f:
            LOG.info("obsm (embeddings): %s", list(f["obsm"].keys()))
        if "raw" in f:
            LOG.info("raw present: yes (raw counts likely in raw/X)")

    # ---- 2) metadata via backed AnnData (loads obs/var/obsm, NOT X/layers) ----
    adata = ad.read_h5ad(args.h5ad, backed="r")
    LOG.info("shape: %s nuclei x %s genes", human(adata.n_obs), human(adata.n_vars))
    LOG.info("X dtype: %s", adata.X.dtype if adata.X is not None else "None")
    LOG.info("obsm keys: %s", list(adata.obsm.keys()))
    LOG.info("obs columns (%d): %s", adata.obs.shape[1], list(adata.obs.columns))

    LOG.info("---- label column summaries ----")
    for col in LABEL_CANDIDATES:
        if col in adata.obs.columns:
            vc = adata.obs[col].value_counts(dropna=False)
            LOG.info("[%s] %d unique; top: %s", col, vc.shape[0],
                     dict(list(vc.head(args.head).items())))

    # quick sanity flags for the integration arm
    emb = [k for k in adata.obsm.keys() if "scanvi" in k.lower() or "scvi" in k.lower()]
    LOG.info("candidate projection embedding(s): %s", emb or "none found — check obsm keys above")

    # ---- 3) optional slim reference for run_heartmap_projection.py ----
    if args.slim:
        import anndata as ad2
        obs = adata.obs.copy()
        var = adata.var.copy()
        obsm = {k: adata.obsm[k][:] for k in adata.obsm.keys()}  # embeddings are small relative to X
        if args.keep_layer:
            # this DOES read a full matrix -> high-mem node required
            LOG.warning("reading layer '%s' into memory (high-mem node needed)", args.keep_layer)
            full = ad.read_h5ad(args.h5ad)  # not backed
            X = full.layers[args.keep_layer] if args.keep_layer in full.layers else full.X
            slim = ad2.AnnData(X=X, obs=obs, var=var, obsm=obsm)
        else:
            import scipy.sparse as sp
            empty = sp.csr_matrix((adata.n_obs, adata.n_vars), dtype="float32")
            slim = ad2.AnnData(X=empty, obs=obs, var=var, obsm=obsm)
        slim.write_h5ad(args.slim)
        LOG.info("slim reference -> %s (%.2f GB)", args.slim, os.path.getsize(args.slim) / 1e9)
        LOG.info("point run_heartmap_projection.py --heartmap_atlas at the dir containing this file.")


if __name__ == "__main__":
    main()
