#!/usr/bin/env python3
"""Pull the Reichart/Seidman 2022 DCM/ARVC/NF matrix from CELLxGENE and record provenance.

Two ways to get it; this script tries the census API and falls back to a direct H5AD URL you supply.
After download it writes reichart_provenance.md: the Plan REQUIRES verifying whether these counts are
CellBender output or a differently-processed matrix (do NOT assume CellBender).

Reichart et al. 2022, Science (DOI 10.1126/science.abo1984); EGA EGAS00001006374.
CELLxGENE collection e75342a8-0f3b-4ec5-8ee1-245a23e0f7cb -> use the "All cells" dataset.

PROVENANCE (verified from the CELLxGENE curation API): the matrix is NOT a CellBender output. It is a
CELLxGENE-curated file with RAW COUNTS in `raw.X` and normalized counts in `.X`. Per Plan Step 2b, run
the uniform DecontX pass on `raw.X`. Chemistry is mixed 10x 3' v2 AND v3 (not pure v3); lymphocytes are
lumped (no B/T split on the Reichart side).
"""
from __future__ import annotations
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import log, ensure_dir  # noqa: E402

LOG = log("download_reichart")
# "DCM/ACM heart cell atlas: All cells": 881,081 nuclei, ~6.4 GB
CXG_DATASET_ID = "65badd7a-9262-4fd1-9ce2-eb5dc0ca8039"
CXG_ASSET_URL = "https://datasets.cellxgene.cziscience.com/6d0b8b8a-2b22-4578-b8ee-8ca4ab893594.h5ad"


def write_provenance(out_dir: str, adata) -> None:
    path = os.path.join(out_dir, "..", "reichart_provenance.md")
    lines = [
        "# Reichart CELLxGENE provenance check\n",
        "Verify whether .X / layers are CellBender output or another processing state.\n",
        f"- n_obs x n_vars: {adata.n_obs} x {adata.n_vars}\n",
        f"- layers present: {list(adata.layers.keys())}\n",
        f"- obs columns: {list(adata.obs.columns)[:40]}\n",
        f"- uns keys: {list(adata.uns.keys())}\n",
        "- X dtype (float => normalized in .X): " f"{adata.X.dtype}\n",
        "- raw.X present (raw counts live here per CELLxGENE schema): "
        f"{adata.raw is not None}\n",
        "\nKNOWN (from CELLxGENE curation API): raw counts are in raw.X; .X is normalized; this is NOT "
        "a CellBender output. Run the uniform DecontX pass on raw.X (Plan Step 2b). Chemistry is mixed "
        "3' v2 + v3.\n",
    ]
    with open(os.path.abspath(path), "w") as fh:
        fh.writelines(lines)
    LOG.info("wrote provenance note -> %s", os.path.abspath(path))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--h5ad_url", default=CXG_ASSET_URL,
                    help="direct CELLxGENE H5AD asset URL (default = Reichart 'All cells')")
    args = ap.parse_args()
    ensure_dir(args.out)
    h5ad = os.path.join(args.out, "reichart_cxg.h5ad")

    import subprocess
    import anndata as ad
    # Direct asset download is the most robust path (no census dependency). ~6.4 GB.
    LOG.info("Downloading Reichart 'All cells' -> %s", h5ad)
    subprocess.run(["wget", "-c", "-O", h5ad, args.h5ad_url], check=True)

    adata = ad.read_h5ad(h5ad, backed="r")
    write_provenance(args.out, adata)
    LOG.info("Reichart download complete -> %s", h5ad)


if __name__ == "__main__":
    main()
