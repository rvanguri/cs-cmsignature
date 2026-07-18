#!/usr/bin/env python3
"""Build the sample manifest every Slurm array job indexes into.

Scans $RAW for raw_feature_bc_matrix files (one per sample), and for each computes an `expected_cells`
SEED from the barcode-rank knee (same algorithm everywhere -> no hand-tuning, per Plan Step 2a).

Output TSV columns (the order 02_cellbender.sbatch relies on):
  1 sample_id  2 dataset  3 accession  4 chemistry  5 procurement  6 region  7 raw_path  8 expected_cells_seed

Disease/procurement/region per sample come from an optional sample_meta.tsv you provide
(sample_id<TAB>dataset<TAB>accession<TAB>chemistry<TAB>procurement<TAB>region). If absent, those cols
are filled "NA" and you should complete them before DE (the DE model needs disease/anatomy/sex).
"""
from __future__ import annotations
import argparse
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import log, require  # noqa: E402

LOG = log("build_manifest")

# filename patterns that indicate a per-sample RAW (unfiltered) droplet matrix
RAW_PATTERNS = ["*raw_feature_bc_matrix.h5", "*raw_feature_bc_matrix",
                "*raw_gene_bc_matrices*", "*unfiltered*.h5"]


def find_raw_matrices(raw_root: str) -> list[str]:
    hits: list[str] = []
    for pat in RAW_PATTERNS:
        hits += glob.glob(os.path.join(raw_root, "**", pat), recursive=True)
    return sorted(set(hits))


# A matrix with at least this many barcodes is treated as RAW (has empty droplets) -> CellBender.
# Fewer barcodes => already cell-called (filtered) -> skip CellBender, DecontX-only (like Reichart).
CELLBENDER_MIN_BARCODES = 20000


def matrix_stats(path: str) -> tuple[int, int]:
    """Return (expected_cells_knee, n_barcodes) for a 10x matrix.

    expected_cells = barcode-rank knee (CellBender --expected-cells seed).
    n_barcodes     = total barcodes in the file (raw matrices have hundreds of thousands; a filtered
                     matrix has only the called cells). Used to decide CellBender vs DecontX-only.
    Robust fallback if the matrix can't be read: (5000, 0).
    """
    try:
        import scanpy as sc
        if path.endswith(".h5"):
            adata = sc.read_10x_h5(path)
        else:
            adata = sc.read_10x_mtx(path if os.path.isdir(path) else os.path.dirname(path))
        n_barcodes = int(adata.n_obs)
        totals = np.asarray(adata.X.sum(axis=1)).ravel()
        totals = np.sort(totals)[::-1]
        totals = totals[totals > 0]
        if totals.size < 50:
            return max(int(totals.size), 100), n_barcodes
        x = np.log10(np.arange(1, totals.size + 1))
        y = np.log10(totals)
        dy = np.gradient(y, x)
        d2y = np.gradient(dy, x)
        knee_idx = int(np.argmin(d2y[: int(0.9 * len(d2y))]))  # avoid tail noise
        return max(knee_idx + 1, 200), n_barcodes
    except Exception as e:  # noqa: BLE001
        LOG.warning("matrix stats failed for %s (%s) -> default (5000, 0)", path, e)
        return 5000, 0


def load_meta(raw_root: str) -> dict[str, list[str]]:
    """Optional sample_meta.tsv keyed by sample_id."""
    meta_path = os.path.join(raw_root, "sample_meta.tsv")
    meta: dict[str, list[str]] = {}
    if os.path.exists(meta_path):
        with open(meta_path) as fh:
            for ln in fh:
                if not ln.strip() or ln.startswith("#"):
                    continue
                p = ln.rstrip("\n").split("\t")
                # sample_id, dataset, accession, chemistry, procurement, region
                meta[p[0]] = (p + ["NA"] * 6)[1:6]
    else:
        LOG.warning("no sample_meta.tsv in %s: dataset/procurement/region will be NA", raw_root)
    return meta


def infer_sample_id(path: str, raw_root: str) -> str:
    rel = os.path.relpath(path, raw_root)
    # use the first path component (dataset/sample/...) or the filename stem
    parts = rel.split(os.sep)
    return parts[1] if len(parts) > 1 else os.path.splitext(parts[0])[0]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    require(args.raw, "raw dir")

    matrices = find_raw_matrices(args.raw)
    if not matrices:
        sys.exit(f"ERROR: no raw_feature_bc_matrix found under {args.raw}")
    LOG.info("found %d raw matrices", len(matrices))
    meta = load_meta(args.raw)

    n_cb = 0
    with open(args.out, "w") as out:
        out.write("sample_id\tdataset\taccession\tchemistry\tprocurement\tregion\traw_path\t"
                  "expected_cells_seed\tn_barcodes\tcellbender\ttotal_droplets\n")
        for m in matrices:
            sid = infer_sample_id(m, args.raw)
            dataset, accession, chem, proc, region = meta.get(sid, ["NA"] * 5)
            knee, n_bc = matrix_stats(m)
            cellbender = 1 if n_bc >= CELLBENDER_MIN_BARCODES else 0
            n_cb += cellbender
            if cellbender:
                # CellBender needs  total_droplets > expected_cells, both <= n_barcodes. The raw knee is
                # unreliable on big curves (can read 100k+), so clamp expected to a sane bound AND keep
                # it well under n_barcodes; total_droplets = expected + buffer guarantees the ordering.
                exp = max(200, min(knee, 20000, n_bc // 3))
                total_droplets = min(n_bc - 1, exp + 25000)
            else:
                exp, total_droplets = knee, 0          # filtered: not used by CellBender
            out.write(f"{sid}\t{dataset}\t{accession}\t{chem}\t{proc}\t{region}\t{m}\t"
                      f"{exp}\t{n_bc}\t{cellbender}\t{total_droplets}\n")
            LOG.info("%s  expected=%d  total_droplets=%d  barcodes=%d  %s", sid, exp, total_droplets, n_bc,
                     "CellBender" if cellbender else "DecontX-only (filtered)")
    LOG.info("manifest -> %s (%d samples; %d CellBender, %d DecontX-only)",
             args.out, len(matrices), n_cb, len(matrices) - n_cb)


if __name__ == "__main__":
    main()
