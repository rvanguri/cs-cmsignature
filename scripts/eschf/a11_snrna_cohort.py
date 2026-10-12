#!/usr/bin/env python3
"""ESC HF figure v3, A11: snRNA CS cohort accounting. Run from repo root.

Joins the Neyazi supplement (Table S1 'Sarcoid' sheet: Tissue Source; Table S2: snRNA samples by anatomical
region LVFW/RVFW/SP and by sarcoid histologic class Preserved/Granuloma/Fibrotic) to the analysed snRNA
cohort (data/snrna/snrna_patient_cohort.tsv: anatomy3, n_specimens, n_cells).
Reports, per S1 patient: whether analysed, Tissue Source, S2 region and histology counts, the analysed
specimen count, and whether S2's two breakdowns agree. Table S2 has no apex category; the sampled region is
reported as S2 records it. The reason a patient is missing is filled from the raw-object check
(results/eschf/a11_raw_donor_check.tsv, BigPurple) when that file exists.
"""
import os
import pandas as pd

s1 = pd.read_excel("data/external/neyazi_table_s1.xlsx", sheet_name="Sarcoid")
s1["Subject"] = s1.Subject.astype(str)
raw = pd.read_excel("data/external/neyazi_table_s2.xlsx", header=None)
hdr = ["Subject", "Diagnosis", "sn_LVFW", "sn_RVFW", "sn_SP", "sn_Preserved", "sn_Granuloma", "sn_Fibrotic",
       "st_LVFW", "st_RVFW", "st_SP", "WGS"]
s2 = raw.iloc[3:, :12].copy(); s2.columns = hdr
s2 = s2[s2.Diagnosis.astype(str).str.contains("Sarcoid")].copy(); s2["Subject"] = s2.Subject.astype(str)
for c in hdr[2:]:
    s2[c] = pd.to_numeric(s2[c])
assert set(s2.Subject) == set(s1.Subject), set(s2.Subject) ^ set(s1.Subject)
coh = pd.read_csv("data/snrna/snrna_patient_cohort.tsv", sep="\t"); coh = coh[coh.disease == "CS"].copy()
coh["Subject"] = coh.patient.str.replace("neyazi_", "", regex=False)
A = s1[["Subject", "Center", "Tissue Source", "LVEF (%)"]].merge(s2.drop(columns="Diagnosis"), on="Subject")
A = A.merge(coh[["Subject", "patient", "anatomy3", "n_specimens", "n_cells"]], on="Subject", how="left")
A["analysed"] = A.patient.notna()
A["s2_n_by_region"] = A[["sn_LVFW", "sn_RVFW", "sn_SP"]].sum(1)
A["s2_n_by_histology"] = A[["sn_Preserved", "sn_Granuloma", "sn_Fibrotic"]].sum(1)
A["s2_breakdowns_agree"] = A.s2_n_by_region == A.s2_n_by_histology
A["specimens_match_s2_region"] = A.n_specimens == A.s2_n_by_region
A["s2_regions_sampled"] = A.apply(lambda r: "+".join(k for k in ("LVFW", "RVFW", "SP") if r[f"sn_{k}"] > 0) or "none recorded", axis=1)
A["apex_recorded"] = "no (Table S2 has no apex category)"
A["missing_reason"] = ""
if os.path.exists("results/eschf/a11_raw_donor_check.tsv"):
    rc = pd.read_csv("results/eschf/a11_raw_donor_check.tsv", sep="\t").set_index("Subject")
    for i, r in A[~A.analysed].iterrows():
        if r.Subject in rc.index:
            A.loc[i, "missing_reason"] = rc.loc[r.Subject, "reason"]
A.loc[(~A.analysed) & (A.missing_reason == ""), "missing_reason"] = "not determined from the supplied tables"
A = A.sort_values(["analysed", "Tissue Source", "Subject"])
A.to_csv("results/eschf/a11_snrna_cohort_accounting.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
print(A[["Subject", "Tissue Source", "analysed", "anatomy3", "s2_regions_sampled", "n_specimens", "s2_n_by_region",
         "s2_n_by_histology", "s2_breakdowns_agree", "specimens_match_s2_region", "missing_reason"]].to_string(index=False))
print(A.groupby(["analysed", "Tissue Source"]).size())
