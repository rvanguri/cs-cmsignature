# A1 — Sections to patients, and patient-level refit (gate)

**Gate**
- Unique CS patients retained: **5** (7 retained sections: CS_100; CS_101-2; CS_102-1, CS_102-2; CS_107; CS_108-1, CS_108-2). 12 CS sections from 9 patients in total (GSE314910).
- Unique comparator patients: **6** (6 sections, six distinct patients, confirmed by the PI).
- The down-in-CS direction survives at patient level. Patient pseudobulk (5 vs 6): 130 genes at FDR 5% (104 down, 26 up). All 192 original hits keep their sign (139/139 down stay down), 127/192 remain FDR < 5%, and 190/192 have a 95% CI excluding 0. Section level with patient blocking (consensus correlation 0.479): 186 hits (136 down, 50 up), 192/192 sign-concordant, 181/192 remain FDR < 5%.
- Patient pseudobulk loses 62 of the 192 hits (32%), which is fewer than "most". Hit counts depend on the unit of analysis, so effect sizes with 95% intervals are the more stable quantity to report alongside them.

**Specification**
- Gene universe fixed to the 985 tested genes. The original section-level fit is reproduced exactly from `data/spatial_first/pseudobulk_T2000.tsv.gz` (asserted in the script).
- Patient pseudobulk: counts summed over retained sections within patient, TMM, voom, `~0 + arm`.
- Patient blocking: voom → `duplicateCorrelation(block = patient)` → voom with block and correlation → second `duplicateCorrelation` → `lmFit(block, correlation)`. Consensus correlation 0.4729 (pass 1) and 0.4790 (pass 2).
- All designs were checked for full rank (rank 2 of 2).
- QC-passing rows restrict the 985-gene fit to the 732 QC-passing genes. FDR is not re-adjusted. Genes outside the locus-QC universe count as not passing.
- 95% CIs use the t quantile on `df.total`. The original table used 1.96, so its intervals are marginally narrower.

| refit | gene set | hits FDR5 (up/down) | r vs original | sign concordance, original hits | original hits still FDR5 |
|---|---|---|---|---|---|
| original section level | all 985 | 192 (53/139) | 1 | 192/192 | 192/192 |
| original section level | QC-pass 732 | 128 (40/88) | 1 | 128/128 | 128/128 |
| patient pseudobulk | all 985 | 130 (26/104) | 0.994 | 192/192 | 127/192 |
| patient pseudobulk | QC-pass 732 | 87 (20/67) | 0.993 | 128/128 | 84/128 |
| section, patient-blocked | all 985 | 186 (50/136) | 0.999 | 192/192 | 181/192 |
| section, patient-blocked | QC-pass 732 | 126 (38/88) | 0.999 | 128/128 | 122/128 |

Values are emitted by `scripts/manuscript_results.py` from `a1_patient_refit_summary.tsv`.

**Granuloma column:** under the A2 definition (immune-dominant spots by lineage argmax), all 12 CS sections contain granulomatous spots (94 or more each). None of the sections lacks a granuloma.
