#!/usr/bin/env python3
"""The two gates (Plan Step 6) + integration validation.

6a integration validation: UMAP by study/disease/cell type; kBET/LISI (scib) if available.
6b DECONTAMINATION QC gate: median log2FC of the non-CM lineage ambient panel inside CS-CMs, per
   contrast. PASS = markers ~0 in every contrast; CS-vs-DCM drops from pre-decontam +1.86 toward ~0.
   Dual-count concordance hook: compare decontaminated vs raw (CellBender vs CellRanger): emits a
   per-gene concordance table (requires both count matrices; TODO wire raw counts path).
6c PROCUREMENT-robustness gate: each candidate must hold across explant-matched mimics
   (DCM, ICM-Simonson, ARVC, Chaffin-LV), not just mismatched (HCM-myectomy, NF).
Writes gate tables (TSV) + figures.
"""
from __future__ import annotations
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import log, require, ensure_dir, PANELS  # noqa: E402

LOG = log("run_gates")


def panel_log2fc(adata, genes, group_mask, ref_mask):
    """median log2FC of a gene panel between two cell groups (pseudo-expression)."""
    genes = [g for g in genes if g in adata.var_names]
    if not genes:
        return np.nan
    import scanpy as sc  # noqa: F401
    sub = adata[:, genes]
    X = sub.X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    g = np.log1p(X[group_mask].mean(axis=0))
    r = np.log1p(X[ref_mask].mean(axis=0))
    return float(np.median(g - r))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--integrated", required=True)
    ap.add_argument("--cellbender_dir", required=True, help="for dual-count concordance (raw counts)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--figdir", required=True)
    ap.add_argument("--disease_col", default="disease")
    ap.add_argument("--celltype_col", default="predicted_cell_type")
    args = ap.parse_args()
    require(args.integrated, "integrated h5ad")
    ensure_dir(args.out); ensure_dir(args.figdir)

    import scanpy as sc
    adata = sc.read_h5ad(args.integrated)

    # ---- 6a: UMAP + integration metrics (visualization only; run on a SUBSAMPLE) ----
    # neighbors/UMAP/scib on ~1M cells take hours and OOM; a 1M-cell UMAP figure is unreadable anyway.
    # The actual gates (6b/6c) below run on the FULL data.
    if "X_scANVI" in adata.obsm:
        import numpy as np
        sub = adata
        if adata.n_obs > 150_000:
            idx = np.random.RandomState(0).choice(adata.n_obs, 150_000, replace=False)
            sub = adata[idx].copy()
            LOG.info("6a: subsampling %d/%d cells for UMAP/metrics", sub.n_obs, adata.n_obs)
        sc.pp.neighbors(sub, use_rep="X_scANVI")
        sc.tl.umap(sub)
        for color in (args.disease_col, "study", args.celltype_col):
            if color in sub.obs:
                sc.pl.umap(sub, color=color, show=False, save=f"_{color}.png")
        try:
            import scib
            metrics = scib.metrics.metrics(sub, sub, "study", args.celltype_col,
                                           embed="X_scANVI", ilisi_=True, kBET_=True)
            metrics.to_csv(os.path.join(args.out, "integration_metrics.csv"))
        except Exception as e:  # noqa: BLE001
            LOG.warning("scib metrics skipped (%s)", e)
        del sub

    # ---- 6b: decontamination QC gate on CS cardiomyocytes ----
    cm = adata.obs[args.celltype_col].astype(str).str.contains("ardiomyocyte", na=False)
    is_cs = adata.obs[args.disease_col].astype(str).str.upper().eq("CS")
    rows = []
    for contrast_disease in ["DCM", "ICM", "ARVC", "HCM", "NF"]:
        ref = adata.obs[args.disease_col].astype(str).str.upper().eq(contrast_disease)
        cs_cm = (cm & is_cs).values
        ref_cm = (cm & ref).values
        if cs_cm.sum() < 10 or ref_cm.sum() < 10:
            continue
        for panel, genes in PANELS["lineage_ambient"].items():
            rows.append({
                "contrast": f"CS_vs_{contrast_disease}", "panel": panel,
                "median_log2fc": panel_log2fc(adata, genes, cs_cm, ref_cm),
            })
    qc = pd.DataFrame(rows)
    qc.to_csv(os.path.join(args.out, "decontam_qc_gate.tsv"), sep="\t", index=False)
    LOG.info("6b decontamination QC gate -> decontam_qc_gate.tsv (%d rows)", len(qc))

    # dual-count concordance: TODO wire raw CellRanger counts and re-test each reported DEG in both;
    # require directional concordance + background heuristic <0.5 (Plan 6b).
    pd.DataFrame({"note": ["wire raw CellRanger counts from --cellbender_dir to complete dual-count "
                           "concordance per Plan Step 6b"]}).to_csv(
        os.path.join(args.out, "dual_count_concordance.TODO.tsv"), sep="\t", index=False)

    # ---- 6c: procurement-robustness scaffold (filled once DE candidates exist) ----
    pd.DataFrame({
        "matched_contrasts": ["CS_vs_DCM;CS_vs_ICM_Simonson;CS_vs_ARVC;CS_vs_HCM_ChaffinLV"],
        "mismatched_contrasts": ["CS_vs_HCM_myectomy;CS_vs_NF"],
        "expected_survive": [";".join(PANELS["cs_cm_expected_up"])],
        "expected_fail_control": [";".join(PANELS["cm_dediff_control"])],
    }).to_csv(os.path.join(args.out, "procurement_gate_spec.tsv"), sep="\t", index=False)

    LOG.info("gates complete -> %s", args.out)


if __name__ == "__main__":
    main()
