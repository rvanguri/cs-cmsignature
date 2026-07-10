#!/usr/bin/env python3
"""Generate a tiny synthetic dataset to dry-run the whole pipeline before real data lands.

Produces (~450 nuclei integrated + ~120 for Liu):
  <out>/raw/<dataset>/<sample>/raw_feature_bc_matrix/{matrix.mtx.gz,features.tsv.gz,barcodes.tsv.gz}
        + real cells over a sea of background droplets (so build_manifest's knee has something to find)
  <out>/raw/sample_meta.tsv                  (disease/procurement/region per sample)
  <out>/decontx/<dataset>_decontx_counts.h5ad  (so QC->scANVI->DE run without the GPU/R legs)
  <out>/raw/liu/liu_qc.h5ad                   (for run_liu_standalone.py)

Signal is wired so the chain finds something: CS cardiomyocytes get boosted TNNI3K/GJB7 (the corrected
expected CS-CM up-genes), CM markers are high in CMs, lineage markers high in their cell types.
NOT biology — just enough structure for every step to execute and produce non-empty outputs.
"""
from __future__ import annotations
import argparse
import gzip
import os
import sys

import numpy as np
import scipy.io
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import log, ensure_dir, CELLTYPE_MARKERS, PANELS  # noqa: E402

LOG = log("make_synthetic")
RNG = np.random.default_rng(0)

# integrated 3' arm: (dataset, sample, disease, anatomy, sex, n_cells)
SAMPLES = [
    ("neyazi",  "cs1",  "CS",   "LV",  "M", 30), ("neyazi", "cs2", "CS", "RV", "F", 28),
    ("neyazi",  "cs3",  "CS",   "LV",  "M", 30),
    ("reichart","dcm1", "DCM",  "LV",  "M", 28), ("reichart","dcm2","DCM","LV","F", 26),
    ("reichart","arvc1","ARVC", "RV",  "M", 26), ("reichart","arvc2","ARVC","RV","F", 24),
    ("reichart","nf1",  "NF",   "LV",  "M", 26), ("reichart","nf2", "NF", "LV", "F", 26),
    ("larson",  "hcm1", "HCM",  "IVS", "M", 26), ("larson", "hcm2","HCM","IVS","F", 24),
    ("chin2022","hcm3", "HCM",  "IVS", "M", 24), ("chin2022","nf3","NF","IVS","F", 24),
    ("chin2021","nf4",  "NF",   "IVS", "M", 24),
]
# Liu 5' standalone arm
LIU_SAMPLES = [
    ("liu", "lcs1",  "CS",  "LV", "M", 30), ("liu", "lcs2", "CS",  "LV", "F", 28),
    ("liu", "licm1", "ICM", "LV", "M", 30), ("liu", "licm2","ICM", "LV", "F", 28),
]
CELL_TYPES = list(CELLTYPE_MARKERS.keys())   # Cardiomyocyte first, then non-CM lineages


def gene_panel() -> list[str]:
    genes: list[str] = []
    for gs in CELLTYPE_MARKERS.values():
        genes += gs
    for key in ("cs_cm_expected_up", "cs_cm_not_expected", "cm_dediff_control",
                "ieg_hsp", "generic_hf_up", "generic_hf_down"):
        genes += PANELS[key]
    genes += PANELS["fib_niche"]["COL22A1_DCM"] + PANELS["fib_niche"]["TNC_ICM"]
    # pad with filler genes so QC/gene-filtering has a realistic width
    genes += [f"FILLER{i}" for i in range(60)]
    genes += ["MT-CO1", "MT-ND1", "MT-CYB"]   # a few MT genes for pct_counts_mt
    # de-dup, preserve order
    seen, out = set(), []
    for g in genes:
        if g not in seen:
            seen.add(g); out.append(g)
    return out


GENES = gene_panel()
GIDX = {g: i for i, g in enumerate(GENES)}


def simulate_cell(cell_type: str, disease: str) -> np.ndarray:
    """Return an integer count vector over GENES for one nucleus."""
    x = RNG.poisson(0.4, size=len(GENES)).astype(int)          # background
    for g in CELLTYPE_MARKERS[cell_type]:                       # this cell type's markers up
        x[GIDX[g]] += RNG.poisson(8)
    if cell_type == "Cardiomyocyte":
        for g in PANELS["cm_dediff_control"]:                  # CM identity genes
            x[GIDX[g]] += RNG.poisson(6)
        if disease == "CS":                                    # CS-CM boosted up-genes (signal)
            for g in PANELS["cs_cm_expected_up"]:
                x[GIDX[g]] += RNG.poisson(10)
        if disease in ("DCM", "ICM", "ARVC", "HCM"):           # failing axis
            for g in PANELS["generic_hf_up"]:
                x[GIDX[g]] += RNG.poisson(5)
    if cell_type == "Fibroblast":
        niche = "COL22A1_DCM" if disease == "DCM" else "TNC_ICM"
        for g in PANELS["fib_niche"][niche]:
            x[GIDX[g]] += RNG.poisson(7)
    return x


def build_anndata(samples):
    import anndata as ad
    import pandas as pd
    rows, obs = [], []
    for dataset, sample, disease, anatomy, sex, n in samples:
        # ~70% cardiomyocytes, rest spread over non-CM lineages
        for k in range(n):
            ct = "Cardiomyocyte" if k < int(0.7 * n) else CELL_TYPES[1 + (k % (len(CELL_TYPES) - 1))]
            rows.append(simulate_cell(ct, disease))
            obs.append(dict(study=dataset, disease=disease, individual=sample,
                            anatomy=anatomy, sex=sex, cell_type=ct,
                            barcode=f"{sample}_{k:03d}"))
    X = sp.csr_matrix(np.vstack(rows))
    obs_df = pd.DataFrame(obs)
    obs_df.index = obs_df["barcode"].values        # set index WITHOUT leaving a clashing column
    obs_df = obs_df.drop(columns=["barcode"])
    obs_df.index.name = None
    var_df = pd.DataFrame(index=GENES)
    return ad.AnnData(X=X, obs=obs_df, var=var_df)


def write_10x_mtx(dest: str, X_cells_by_genes: sp.spmatrix, barcodes: list[str]) -> None:
    """Write a 10x v3 raw matrix (features x barcodes) gzipped, readable by scanpy.read_10x_mtx."""
    ensure_dir(dest)
    feats = (X_cells_by_genes.T).tocoo()                        # genes x cells
    tmp = os.path.join(dest, "matrix.mtx")
    scipy.io.mmwrite(tmp, feats, field="integer")
    with open(tmp, "rb") as fi, gzip.open(tmp + ".gz", "wb") as fo:
        fo.writelines(fi)
    os.remove(tmp)
    with gzip.open(os.path.join(dest, "features.tsv.gz"), "wt") as f:
        for g in GENES:
            f.write(f"ENSG_{g}\t{g}\tGene Expression\n")
    with gzip.open(os.path.join(dest, "barcodes.tsv.gz"), "wt") as f:
        for b in barcodes:
            f.write(b + "\n")


def write_raw_with_background(out_root: str) -> None:
    """Per-sample raw droplet matrix = real cells + many low-count background barcodes (knee)."""
    for dataset, sample, disease, anatomy, sex, n in SAMPLES:
        cells = np.vstack([simulate_cell("Cardiomyocyte" if k < int(0.7 * n)
                                         else CELL_TYPES[1 + (k % (len(CELL_TYPES) - 1))], disease)
                           for k in range(n)])
        # real cells get extra depth so they sit above the knee
        cells = cells + RNG.poisson(3, size=cells.shape)
        bg = RNG.poisson(0.15, size=(1500, len(GENES)))         # background droplets
        X = sp.csr_matrix(np.vstack([cells, bg]).astype(int))
        barcodes = [f"{sample}_C{k:03d}" for k in range(n)] + [f"{sample}_B{j:04d}" for j in range(1500)]
        dest = os.path.join(out_root, dataset, sample, "raw_feature_bc_matrix")
        write_10x_mtx(dest, X, barcodes)
    LOG.info("wrote raw 10x mtx for %d samples", len(SAMPLES))


def write_sample_meta(out_root: str) -> None:
    proc = {"CS": "explant", "DCM": "explant", "ICM": "explant", "ARVC": "explant",
            "HCM": "myectomy", "NF": "donor"}
    path = os.path.join(out_root, "sample_meta.tsv")
    with open(path, "w") as f:
        f.write("sample_id\tdataset\taccession\tchemistry\tprocurement\tregion\n")
        for dataset, sample, disease, anatomy, sex, n in SAMPLES:
            chem = "5p" if dataset == "liu" else "3p"
            f.write(f"{sample}\t{dataset}\tSYNTH\t{chem}\t{proc[disease]}\t{anatomy}\n")
    LOG.info("wrote %s", path)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, help="synthetic project root")
    args = ap.parse_args()
    raw = ensure_dir(os.path.join(args.out, "raw"))
    dec = ensure_dir(os.path.join(args.out, "decontx"))
    ensure_dir(os.path.join(raw, "liu"))

    # 1) raw 10x mtx (+ background) for manifest/knee + optional CellBender/DecontX legs
    write_raw_with_background(raw)
    write_sample_meta(raw)

    # 2) decontx-level h5ads per dataset (entry point for the CPU chain: QC->scANVI->DE)
    by_ds: dict[str, list] = {}
    for s in SAMPLES:
        by_ds.setdefault(s[0], []).append(s)
    for ds, samps in by_ds.items():
        adata = build_anndata(samps)
        out = os.path.join(dec, f"{ds}_decontx_counts.h5ad")
        adata.write_h5ad(out)
        LOG.info("%s_decontx_counts.h5ad: %d nuclei x %d genes", ds, adata.n_obs, adata.n_vars)

    # 3) Liu QC'd h5ad (standalone 5' arm) — give it predicted_cell_type for run_liu's CM filter
    liu = build_anndata(LIU_SAMPLES)
    liu.obs["predicted_cell_type"] = liu.obs["cell_type"]
    liu_out = os.path.join(raw, "liu", "liu_qc.h5ad")
    liu.write_h5ad(liu_out)
    LOG.info("liu_qc.h5ad: %d nuclei x %d genes -> %s", liu.n_obs, liu.n_vars, liu_out)

    LOG.info("synthetic data ready under %s", args.out)


if __name__ == "__main__":
    main()
