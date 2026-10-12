# ESC HF v3: discovery pipeline rerun without the Kuppe normal sections

> Record of the rerun, kept as written. Since then the `ESCHF_V3` switch has been removed and this is the only path. `run_v3.sh` and `v3_derived_inputs.py` are now `run_pipeline.sh` and `derived_inputs.py`, and the default-mode control (section 7) no longer applies. See `docs/DISCREPANCIES.md` section 3.

Results are reported as they fall. Interpretation and framing are left to the authors.

## 1. What changed

- The four Kuppe non-failing ("NORMAL_") sections were removed from the discovery export. The export now pools 18 sections (12 CS, 6 comparator) instead of 22. Every other specification is unchanged: 2,000-UMI thinning, five lineage gene sets, argmax assignment, 50-spot floor, 25 % detection rule, limma-voom design, the ±3 snRNA screen, null procedures and seeds.
- Consequences the change forces (agreed before running):
  - **Gene space.** The inner join is now taken over the 18 remaining sections: 17,619 genes (committed) → 18,039 genes (v3).
  - **Q2 dropped.** Q2 (CS / comparators vs normal) and the ‡ category require the normal sections, so neither is computed.
  - **Locus QC.** Same five flags, now from two sources (foong, new_visium). Each fit's own tested genes form the universe (`03_locus_qc.py --universe`).
- Every code change sits behind `ESCHF_V3=1`. Without the variable the edited scripts reproduce the committed results (section 7).

## 2. How it was run

| round | job (BigPurple) | content |
|---|---|---|
| 1 | 6483bd0c-bcb0-4b4d-835d-fc703bcd96d3 | v3 discovery export (`00_export_from_raw.py`, ESCHF_V3=1) |
| 1 | 76a74f35-d5f5-4218-90e2-590c73b01a7f | control: same export with Kuppe, ESCHF_V3 unset; reproduced the committed pseudobulk/region tables |
| 1 | cd11d3d3-433b-4cf0-af46-7f06cf92ff5c | diagnostic: Kuppe removed from the z-pool only, old gene space kept |
| 2 | e87cad27-382f-4ca2-b869-ce481e9b08a5 | primary all-cardiomyocyte arm export (`a0_primary_export.py`) |
| 3 | 6815f9ae-dd83-4a00-8aae-12fdaea62e62 | per-spot program-score export (`a2a6_spot_export.py`) |

The planned two cluster rounds became three. The per-spot export needs the v3 primary program, which comes from the v3 sensitivity fit, so it could not run until that fit existed. The local cascade ran in an isolated copy of the repo. `run_v3.sh` replays the whole sequence (stages cluster1 … compare).

Output layout: the regenerated files sit under `results/eschf/primary_v3/` with the repo's own relative paths (`data/`, `metadata/`, `results/`, `figures/`). The v3 primary-arm files keep the existing `primary_v2` path names, because `eschf_paths.py` only knows the `old` and `v2` arm modes. No file or column was renamed. For the same reason the tested-gene list keeps its committed file name, `metadata/tested_genes_985.txt`, even though in v3 it holds 990 genes.

## 3. Spot and section universe (T = 2,000)

- Spots passing thinning: 61,569 → 61,645.
- Spots whose argmax lineage changed: 29.3 %.
  - Kuppe removed from the z-pool alone (old gene space; the diagnostic job): 1.75 %.
  - Re-thinning over the new gene space (diagnostic → v3): 29.5 %.
- Primary all-CM arm: no section crossed the 50-spot floor.
- ≥ 0.5 mm sensitivity arm:
  - CS_100 (61 → 27 spots) and CS_102-2 (80 → 35) fell below the floor.
  - CS_106 (21 → 60) rose above it.
  - Evaluable CS sections: 7 → 6.
  - Patient 108 is now the only CS patient with more than one evaluable section.

Tables: `lineage_change_by_section_T2000.tsv`, `lineage_change_decomposition.tsv`, `lineage_transition_T2000_old_vs_v3.tsv`, `lineage_change_report.png`.

## 4. Headline values, committed → v3

Sensitivity arm (≥ 0.5 mm; spatial-first manuscript):

| quantity | committed | v3 |
|---|---|---|
| genes tested at 25 % detectability | 985 | 990 |
| FDR < 5 % hits (lower / higher in CS) | 192 (139 / 53) | 193 (133 / 60) |
| hits passing locus QC | 128 | 153 |
| hits stable in every leave-one-section-out refit | 141 | 123 |
| hits stable in every leave-one-disease-out refit | 117 | 105 |
| LOO-stable lower core (file `q1_stable_down_<n>.tsv`) | 115 | 105 |
| core genes testable in snRNA | 134 | 116 |
| concordant and snRNA FDR < 5 %, lower / higher in CS | 15 / 3 | 12 / 1 |
| sign concordance among lower-expression genes | 60 % | 57 % |
| effect-magnitude correlation across platforms | rho = 0.03 (p = 0.76) | rho = −0.06 (p = 0.52) |
| published comparison, DCM: lower in disease; binomial p | 25/36; 0.009 | 22/32; 0.016 |
| published comparison, HCM: lower in disease; binomial p | 22/35; 0.023 | 20/33; 0.048 |
| named metabolic/ion/adhesion genes all lower in CS | all 9 | SLC2A4 no longer in the set |
| MT-free renormalised control | 192/192 same direction, 190/192 FDR < 5 % | 190/190 same direction, 184/190 FDR < 5 % |
| panel-matched comparators | 176/176 same direction | 180/180 same direction |
| QC-pass CS-lower program (`metadata/geneset_primary_program_qcpass_down.txt`) | 88 genes | 101 genes |

Primary all-CM arm (ESC HF revision):

| quantity | committed | v3 |
|---|---|---|
| tested genes | 982 | 985 |
| FDR < 5 % hits (lower in CS) | 222 (150) | 223 (151) |
| QC-passing hits | 148 | 180 |
| program genes (QC-passing, lower in CS) | 95 | 116 |
| hits shared with the sensitivity fit | 162 | 158 |
| A2 spot-level distance CI excludes 0: program / ACADVL / OXPHOS / CM-structural | no / yes / no / yes | no / no / yes / yes |
| A2 distance beta (statsmodels (1\|section)): program / ACADVL | −0.0038 / 0.0982 | −0.0036 / 0.0384 |
| A2 per-patient mean slope: program / ACADVL | 0.0036 / 0.0156 | 0.0209 / 0.0624 |
| A2 patients / sections / spots in model | 8 / 11 / 16,435 | 8 / 11 / 16,125 |
| A4 Chaffin DCM, QC-passing: lower / testable; binomial p | 15/17; 0.0007 | 19/24; 0.002 |
| A4 Chaffin HCM, QC-passing: lower / testable; binomial p | 19/24; 0.0006 | 20/28; 0.0038 |
| A6 immune-dominant fraction, Mann-Whitney p | 0.036 | 0.036 |
| a7 Figure panel D genes | 17 | 19 (adds LMO7, S100A13) |

The A2 and A6 rows have identical values in the sensitivity block of the harness.

## 5. Mechanical change flags (`comparison/old_vs_new_manuscript_values.tsv`)

The v3 run evaluates 878 claims: 850 OK, 23 MISMATCH against the committed manuscript values, 3 NOT_RUN and 2 NOT_REPRODUCIBLE (the same two as in the committed run). Taking the union with the committed claims gives 881 rows:
- 870 in both runs, 9 only in v3, 2 only in the committed run;
- 587 numeric values changed;
- 27 status changes;
- 7 sign changes;
- 0 claims whose p / FDR value crosses 0.05. The closest is the published-comparison HCM binomial p, 0.023 → 0.048.

Sign changes (numeric claims):
- a14 fatty-acid degradation mean slope / mm: 0.042 → −0.0465
- a14 TCA cycle mean slope / mm: −0.0176 → 0.0539
- a14 TCA cycle largest leave-one-gene-out change: −0.045 → 0.0423
- a14 OXPHOS null 2.5 % quantile: 0.0022 → −0.058
- A3 down, QC-passing only, Spearman rho: −0.0099 → 0.0346
- Fig C screened-down rho: −0.0096 → 0.0318
- A5b LVEF Spearman rho, program score: 0.0113 → −0.0068

Categorical changes beyond the A2 CI rows above:
- A10 TCA cycle, 1–2 mm bin: the CS mean is now below the comparator minimum (False → True).
- A9 3a failures in the section-stable core panel: NDUFC2 → ENG.
- The Fig 1C triangle list is unchanged: five screened genes (HSPB2, CALM2, NDUFC2, ENG, CFD). CHMP3 (+2.9) is plotted in range on the widened axis (see section 6).
- Fig D marks:
  - The double dagger (outside QC universe) changes 1 → 0 for ACSL1, ARHGAP1, FBXW5, FITM2 and JUP.
  - The dagger (locus outlier) changes 1 → 0 for IDH2.
  - SPTAN1 keeps its dagger.
- Fig B per-patient slope CIs:
  - program: −0.0266 to 0.0338 → −0.0211 to 0.0629;
  - ACADVL: −0.2341 to 0.2653 → −0.2118 to 0.3367.
  - Both intervals still include 0.

The 9 v3-only and 2 committed-only claims have two sources:
- The A2 per-patient gap-share claims embed the maximum distance in their names (1.73 mm committed, 1.64 mm v3), so each run's version appears as a separate row.
- New panel-D and triangle genes generate claims of their own.

## 6. Stops, deviations and not-regenerated items

1. **A1 refit 2 (section-level, duplicateCorrelation patient block), sensitivity arm: NOT ESTIMABLE.** Only patient 108 still contributes more than one section, and the consensus correlation is NaN. No substitute model was fitted. a1 writes NA, a8 carries it as a nullable integer, and the harness reports the two claims as NOT_RUN. The primary-arm refit is unaffected (203 / 161 hits).
2. **Cardiomyocyte-purity control: NOT RUN.** The `q1_purityadj` fit needs `cm_purity`, which no script or data file in the repo defines. The claim "result persists under adjustment for cardiomyocyte content" is NOT_RUN.
3. **Figure C panel assert.** a7 asserts that the genes beyond the panel C axis are exactly the screen-3a failures. In v3 one QC-passing primary hit, CHMP3 (pooled snRNA log2FC 2.9), falls between the axis half-range (2.5) and the screen threshold (3). In the first v3 pass the axis was left at ±2.5, the assert was gated, and CHMP3 was drawn as a triangle. At Rami's request the panel C y-axis was then widened to ±3 under v3 only (`YLIM_C = 3.0 if _V3 else 2.5`; the matching values in a8 and the harness follow it), so the axis equals the screen. With the axis at ±3 the assert is unconditional again and passes. Panel C shows 170 points and 5 triangles, exactly the five screen-3a failures, and the harness claim "Fig C off-axis genes = screen-3a failures" is OK. The screen itself (|log2FC| > 3) is unchanged. In default mode the axis stays at ±2.5, and a7, a8 and the figures are byte- and pixel-identical to the committed outputs.
4. **A13 (GPNMB) not rerun.** Its inputs come from the Kuppe-inclusive pipeline. Under v3 the 31 A13 claims carry the note "A13 tables are committed values". The A13 spot-count consistency claim is a MISMATCH because the primary spot counts changed.
5. **Spatial-first Figure 1 panel tables not regenerated.** `figure1_panelA_values.tsv` and `figure1_panelB_values.tsv` have no generator in the repo; `figure1_build.py` only renders them. The three dependent claims are marked as committed values, and `Figure1.png` was not redrawn. Panel A is described as the core genes testable in both arms, which number 134 in the committed run and 116 in v3.
6. **Stale stable-core files.** Files whose names embed a core size (`q1_stable_down_115.tsv` etc.) remain from the committed run. The v3 core file is `q1_stable_down_105.tsv` in the sensitivity arm and `q1_stable_down_128.tsv` in the primary arm. The harness resolves concordance labels by prefix.
7. **Locus QC annotation.** One Ensembl REST call returned HTTP 500 and succeeded on rerun. Retired symbols receive the existing missing-annotation flag.
8. **Missing genes.** The `a2`, `a5` and `a5b` missing-gene lists are identical to the committed ones (1, 3, 1 genes; MYH6 is absent from the per-spot export, as before).
9. **Harness edits under v3:**
   - the Q2 read is skipped;
   - a NOT_RUN status was added;
   - concordance labels that embed counts are looked up by prefix;
   - one A2 singular-fit assert is relaxed to its variance bound, with the flag recorded in the note;
   - stale-input notes were added (items 4–5).

## 7. Verification

- **Default mode (ESCHF_V3 unset) in a separate copy.** The cascade of edited local scripts and the harness reproduce the committed run: 870 claims OK, 0 MISMATCH, 2 NOT_REPRODUCIBLE.
  - Regenerated tables are byte-identical to the committed files, except `a0_vs_old_primary.tsv`, which differs in the last floating-point digit of correlation columns.
  - The figures differ only in gene-label placement; the plotted points are identical.
- **Coverage gaps in the default-mode check.** The unedited 02 LOO, 04, 06 and primary-arm locus QC steps were not rerun in default mode. The edited cluster scripts `a0_primary_export.py`, `a2a6_spot_export.py` and `a13_gpnmb_remote.py` were not run in default mode either. `00_export_from_raw.py` was: its round-1 control job matched.
- **Cluster export checks.** The round-2 and round-3 exports reproduced the v3 discovery region and spot counts (`verify_region_counts.tsv`).

## 8. Files

- `code_changes.diff`: every code change of this session, against the pre-session working copy. That baseline equals git HEAD for all tracked files except `scripts/manuscript_results.py`, whose pre-session edits are kept in the baseline. 20 edited files and 3 new ones (`v3_derived_inputs.py`, `v3_old_vs_new.py`, `run_v3.sh`).
- `run_v3.sh`: full reproduction script (also at `scripts/eschf/run_v3.sh`).
- `comparison/old_vs_new_manuscript_values.tsv`: per-claim committed vs v3 values, with the flags above.
- `comparison/old_vs_new_gene_lists.tsv`: 41 gene lists (30 differ), with shared / only-old / only-v3 members.
- `results/manuscript_results.tsv`, `comparison/harness_v3.log`: v3 harness output.
- Nothing was committed to git, and no committed result was overwritten.
