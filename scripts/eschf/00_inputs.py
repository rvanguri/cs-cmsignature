#!/usr/bin/env python3
"""ESC HF short report, step 0: stage inputs and build the annotated hit table.

Run from the repo root.

1. results/eschf/q1_hits_annotated.tsv: the 192 FDR<5% hits of the primary contrast
   (q1_raw_T2000_d25.tsv) joined to the locus QC table (q2_locus_qc.tsv). Locus QC was
   evaluated on the genes tested in the CS-vs-normal contrast, so some hits fall outside
   its universe; for those qc_pass is NA and the gene is NOT counted as QC-passing (same
   convention as scripts/manuscript_results.py).
   flag_length_deviation is renamed flag_locus_quant_outlier: a locus-quantification
   outlier flag (mean expression regressed on log gene span within each dataset; flagged if
   the residual exceeds a fixed margin in any dataset). It is computed on expression level,
   not on effect size.

2. data/snrna/snrna_cm_patient_counts.tsv.gz: patient-level cardiomyocyte pseudobulk counts
   (specimens summed within patient) for the 121 patients of snrna_patient_cohort.tsv.
   Built from the specimen-level CM pseudobulk (--pseudobulk_dir). Per-patient nucleus and
   specimen counts match the cohort table exactly, and limma-voom ~0 + disease + anatomy3
   (filterByExpr by disease, TMM) on this matrix reproduces data/snrna/snrna_CS_vs_*.tsv to
   machine precision (asserted in scripts/eschf/a5b_procurement.R).

3. data/external/neyazi_clinical.tsv: per-patient CS clinical table (Tissue Source, LVEF)
   from the Neyazi supplement Table S1 ("Sarcoid" sheet).
"""
import argparse
import os

import pandas as pd

RES = "results/spatial_first"
OUT = "results/eschf"
FDR = 0.05


def build_hits():
    q1 = pd.read_csv(f"{RES}/q1_raw_T2000_d25.tsv", sep="\t")
    qc = pd.read_csv(f"{RES}/q2_locus_qc.tsv", sep="\t")
    cols = ["gene", "qc_pass", "flag_length_deviation", "flag_segdup",
            "flag_missing_annotation", "flag_readthrough", "flag_source_specific_absence",
            "max_abs_length_dev_log2", "span_kb", "biotype"]
    h = q1[q1["adj.P.Val"] < FDR].merge(qc[cols], on="gene", how="left")
    h["in_qc_universe"] = h.qc_pass.notna()
    h["qc_pass"] = h.qc_pass.astype("boolean").fillna(False).astype(bool)
    h["direction"] = (h.logFC > 0).map({True: "up", False: "down"})
    h = h.rename(columns={"flag_length_deviation": "flag_locus_quant_outlier"})
    h.sort_values("logFC").to_csv(f"{OUT}/q1_hits_annotated.tsv", sep="\t", index=False)
    return h


def build_patient_counts(pbdir):
    meta = pd.read_csv(f"{pbdir}/cm_pseudobulk_meta.tsv", sep="\t")
    cnt = pd.read_csv(f"{pbdir}/cm_pseudobulk_counts.tsv", sep="\t", index_col=0)
    coh = pd.read_csv("data/snrna/snrna_patient_cohort.tsv", sep="\t")
    meta["patient"] = [("neyazi_" + "_".join(i.split("_")[1:-2])) if s == "neyazi" else f"{s}_{i}"
                       for s, i in zip(meta.study, meta.individual)]
    meta = meta[meta.patient.isin(coh.patient)]
    agg = meta.groupby("patient").agg(n_cells=("n_cells", "sum"), n_spec=("individual", "size"))
    chk = coh.set_index("patient").join(agg, rsuffix="_pb")
    assert (chk.n_cells == chk.n_cells_pb).all() and (chk.n_specimens == chk.n_spec).all()
    pc = cnt[meta.individual].T.groupby(meta.set_index("individual").patient).sum().T[coh.patient]
    pc.to_csv("data/snrna/snrna_cm_patient_counts.tsv.gz", sep="\t")
    return pc.shape


def build_neyazi_clinical():
    s1 = pd.read_excel("data/external/neyazi_table_s1.xlsx", sheet_name="Sarcoid")
    keep = ["Subject", "Center", "Gender", "Age", "Tissue Source", "Immunosuppression", "LVEF (%)"]
    s1 = s1[keep].copy()
    s1["patient"] = "neyazi_" + s1.Subject.astype(str)
    # LVEF recorded as '<15' is left-censored; numeric value kept NA, raw string retained
    s1["lvef_raw"] = s1["LVEF (%)"].astype(str)
    s1["lvef"] = pd.to_numeric(s1["LVEF (%)"], errors="coerce")
    s1.to_csv("data/external/neyazi_clinical.tsv", sep="\t", index=False)
    return s1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pseudobulk_dir", default=None,
                    help="dir with cm_pseudobulk_counts.tsv / cm_pseudobulk_meta.tsv")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    h = build_hits()
    print(f"[hits] {len(h)}: {int((h.direction == 'down').sum())} down / "
          f"{int((h.direction == 'up').sum())} up; QC universe {int(h.in_qc_universe.sum())}; "
          f"QC-pass {int(h.qc_pass.sum())} ({int((h.qc_pass & (h.direction == 'down')).sum())} down)")
    if a.pseudobulk_dir:
        print("[snrna patient counts]", build_patient_counts(a.pseudobulk_dir))
    c = build_neyazi_clinical()
    print("[neyazi clinical]", c["Tissue Source"].value_counts().to_dict())


if __name__ == "__main__":
    main()
