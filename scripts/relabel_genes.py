#!/usr/bin/env python3
"""Relabel a .h5ad's var_names from Ensembl IDs to HGNC symbols, IN PLACE.

The pipeline's decontx step used Ensembl gene IDs as var_names, but the gene panels/markers
use gene symbols: so panel/marker/gene-overlap operations silently matched nothing.
This maps Ensembl -> symbol using a 10x features.tsv.gz (col1=Ensembl, col2=Symbol) and renames, keeping
the old id in var['ensembl']. Lets us fix already-computed outputs without re-running scANVI.

Usage:
  relabel_genes.py --h5ad <file.h5ad> --features <features.tsv.gz> [--out <file.h5ad>]
"""
from __future__ import annotations
import argparse
import gzip
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import log, require  # noqa: E402

LOG = log("relabel_genes")


def load_map(features: str) -> dict:
    op = gzip.open if features.endswith(".gz") else open
    m = {}
    with op(features, "rt") as fh:
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            if len(p) >= 2 and p[0] and p[1]:
                m[p[0].upper()] = p[1]        # Ensembl (upper) -> symbol
    return m


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--features", required=True, help="10x features.tsv.gz (Ensembl<TAB>Symbol<TAB>type)")
    ap.add_argument("--out", default="")
    args = ap.parse_args()
    require(args.h5ad, "h5ad"); require(args.features, "features.tsv")

    import anndata as ad
    import pandas as pd

    m = load_map(args.features)
    LOG.info("ensembl->symbol map: %d entries", len(m))
    adata = ad.read_h5ad(args.h5ad)
    ens = adata.var_names.astype(str)
    adata.var["ensembl"] = ens.values
    sym = pd.Series([m.get(e.upper(), e) for e in ens], index=range(len(ens)))
    # disambiguate duplicate symbols with the ensembl id (plan Step 4: "resolve duplicates w/ suffix")
    dup = sym.duplicated(keep=False)
    if dup.any():
        sym = sym.where(~dup.values, sym.astype(str) + "." + ens.values)
    adata.var_names = sym.values
    adata.var_names_make_unique()

    out = args.out or args.h5ad
    adata.write_h5ad(out)
    n_mapped = sum(1 for e in ens if e.upper() in m)
    LOG.info("relabeled %d/%d genes to symbols -> %s", n_mapped, len(ens), out)


if __name__ == "__main__":
    main()
