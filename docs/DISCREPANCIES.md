# Verification: what reproduces, and what it cannot

`python scripts/manuscript_results.py` checks every number in the ESC HF short report (*Cardiomyocyte
transcription in cardiac sarcoidosis*) against values recomputed from the committed tables. It exits non-zero on any
mismatch. As of this commit **865 claims reproduce, none mismatch, and 2 are NOT_RUN** (section 1).

| Block | Claims | What it checks |
|---|---|---|
| `ShortReport` | 62 | each number in the manuscript text and legend, as printed, recomputed and compared at printed precision |
| `ESCHF[primary]` | 114 | A1–A6 values for the primary fit (all cardiomyocyte-dominant spots) |
| `ESCHF[sens>=0.5mm]` | 114 (2 NOT_RUN) | the same checks for the ≥0.5 mm exclusion fit |
| `ESCHF[rev]` | 217 | figure-revision values (per-patient slopes, per-spot scores, panel tables) |
| `ESCHF[v3]` | 108 | figure v3 values (snRNA screen, A10, A11) |
| `ESCHF[a14]` | 252 | pathway scores with fit hits removed |

`bash scripts/eschf/run_pipeline.sh all` regenerates every table from the committed exports and finishes with
the harness.

## 1. Not run (2)

`ESCHF[sens>=0.5mm] A1 all_985 dupcor-blocked hits` and `... qc_pass dupcor-blocked hits`. These are the
patient-blocked (duplicateCorrelation) refits of the ≥0.5 mm fit. Only one patient contributes more than one section
to that fit, and the consensus within-patient correlation comes out NaN. `a1_patient_refit.R` records the refit as
NOT ESTIMABLE instead of fitting it, and the harness reports NOT_RUN. The primary fit's patient-blocked refit is run
and checked.

## 2. Earlier reports

The spatial-first brief report (*Metabolic gene expression differs in non-lesional myocardium between cardiac
sarcoidosis and other cardiomyopathies*), with its driver `scripts/spatial_first/run_all.sh`,
`scripts/figure1_build.py`, `figures/Figure1.*`, the figure 1 panel tables, its 40-claim harness block and this file's
earlier account of it, is preserved at tag `brief-report-spatial-first`. The superseded single-nucleus arm is at tag
`pre-slim-2026-09-14`.

## 3. This commit: the analysis without normal donors

- **Normal-donor Visium removed.** The four Kuppe et al. non-failing sections are no longer read anywhere. The
  discovery export pools the 12 CS and 6 comparator sections. The Q2 contrast (CS and comparators vs normal), the ‡
  category, the purity-adjusted fit and the GS22 diagnostic required those sections and are removed. Locus QC uses two
  sources, with each fit's tested genes as its universe. Earlier this was the `ESCHF_V3=1` path; the switch is gone
  and this is now the only path. Old-vs-new values, gene lists and the code diff are in `results/eschf/primary_v3/`.
- **A13 removed.** The GPNMB source check depended on the normal sections. Its scripts, tables and its section in
  `results/eschf/SUMMARY.md` are removed.
- **Manuscript-text block.** `short_report_claims()` adds the 62 `ShortReport` claims. The explant-only expected
  agreement for lower-in-CS genes is 49% (`a5b_txonly_concordance.tsv`, new `null_mean` column). An earlier text
  printed 48%, which is the full-cohort background; the text was corrected to 49%.
- **Locus-QC annotation is fully cached.** `03_locus_qc.py` reads `results/spatial_first/q2_locus_qc_ensembl_fetched.tsv`
  (90 genes absent from `data/spatial_first/ensembl_annotation_tested.tsv`) and
  `q2_locus_qc_ensembl_unresolved.txt` (14 symbols Ensembl does not resolve, flagged `flag_missing_annotation`), so a
  rerun needs no network. Before this, the 90 genes were re-queried on every run. `q2_locus_qc.tsv` is byte-identical
  either way.
- **Figure label placement is not deterministic.** Rerunning `a7_figure.py` can move a few gene labels in panel C.
  Plotted values and `a7_panel_values.tsv` do not change.
- **Pipeline driver.** `scripts/eschf/run_pipeline.sh` (formerly `run_v3.sh`) and `scripts/eschf/derived_inputs.py`
  (formerly `v3_derived_inputs.py`).

Sections 4–4.5 below are the revision record kept as written. They predate the rerun without normal donors, so their
counts are pre-rerun values (for example, 17 genes in panel D, now 19). Current values are in
`results/eschf/SUMMARY.md` and `results/manuscript_results.tsv`.

## 4. ESC HF short report (A1–A7)

`eschf_claims()` in `scripts/manuscript_results.py` checks each value that the short report can quote
against an independent recomputation from upstream inputs: the per-spot export, the annotated hit table,
the snRNA contrast and pooled tables, the published Chaffin tables, the Hamed G2−G1 table and the CS
patient-score table. It runs the same A1–A6 checks once per fit:

- `ESCHF[primary]` (114 claims): the primary contrast over **all** CM-dominant spots (figure-fixes handoff
  s.0); tables in `results/eschf/primary_v2/`, inputs in `data/eschf/primary_v2/`.
- `ESCHF[sens>=0.5mm]` (114 claims): the earlier fit on CM-dominant spots ≥0.5 mm from an immune-dominant
  spot, now a sensitivity analysis; tables at the top level of `results/eschf/`. These are the same 114
  checks as before the revision, with identical recomputed values.
- `ESCHF[rev]` (211 claims): the revision items — new vs earlier fit (hit counts, overlap, sign retention,
  program gene list), the per-patient slopes and gap shares (refitted with statsmodels OLS from the per-spot
  export), the OXPHOS level table, the descriptive direction table, and every Figure 1 (v3, four panels) and
  legend number in `results/eschf/a7_panel_values.tsv` (v5 figure): Panel A spot counts; Panel B score-gene count, bin
  means, slopes and max distance; Panel C plotted, off-axis and screen-3a counts, each off-axis triangle's value
  and its open neutral drawing (v4), screened and unscreened sign agreement, ρ and background; Panel D gene set (section-stable core
  recomputed from the six leave-one-comparator-section-out refits) and every †/‡/§ mark, the snRNA vs NF magnitude bound; the
  patient-pseudobulk sentence, the v4 Panel B pathway sentence (against A14) and the v5 Panel C source-absence
  sentence (3b rule recomputed from per-study mean log2 CPM); fitted section and patient counts.
- `ESCHF[v3]` (139 claims): the v3 analyses (`HANDOFF_figure_v3_ESCHF.md`), recomputed from upstream inputs —
  A9 3a/3b gene sets and counts per set, the Chaffin |log2FC| > 3 counts, the concordant sets (all hits, after
  3a, core), the sign-agreement counts; A12 KCTD14 absence; A10 spot-level bin scores from the spot count
  matrix and patient-level differences and Welch p from the voom-form patient log2 CPM, plus member/hit
  overlap; A11 unanalysed subject, specimen-count match against raw Table S2, and the GEO listing; A13
  spatial GPNMB summaries from the per-spot export (with a check that the remote GPNMB counts equal the
  committed spot counts), snRNA subject-level comparisons, per-patient distance correlations, and GPNMB
  passing both screens. Permutation backgrounds are checked against `a9_concordance_screened.tsv`, not
  regenerated.
- `ESCHF[a14]` (252 claims): A14 (`HANDOFF_pathway_hits_removed.md`), recomputed from the hit tables, KEGG
  files, universe, voom-form patient log2 CPM, spot counts and model input — members and removed genes of each
  version, too-small flags, patient-level means, ranges, Welch intervals and p (scipy Welch), the null pool,
  decile composition (explicit quantile edges) and, for every one of the 10,000 draws per pathway, size,
  distinct genes, pool membership and decile match, with each draw's statistic recomputed to give the null
  mean, quantiles, percentile and p; leave-one-gene-out range and most influential gene, per-member
  direction counts, distance-bin means from the spot count matrix, per-patient slopes (statsmodels OLS); and
  the SUMMARY.md Figure 1C triangle list (five genes, values equal to the pooled snRNA log2FC).

All 319 agree (tolerance 1e-6, except the mixed-model slopes: 1%).

**What is recomputed.** A1 retained sections and patients per arm, hit counts and sign retention for the
section-level, patient-pseudobulk and duplicateCorrelation fits. A2 spot and patient counts per distance bin,
the mean of patient means per bin for the program score and *ACADVL*, the maximum CM-spot distance, the
comparator reference (mean, min, max of patient means), and the distance slopes for the program score,
*ACADVL* and the positive control. A3 testable counts, sign-agreement fractions and Spearman ρ, all hits vs
QC-pass, both directions. A4 Chaffin testable and lower-in-disease counts and binomial p; internal CS-vs-NF
testable count and fraction lower. A5 overlap size, the severity-axis Spearman ρ, and the G2-lower/G2-higher
counts among CS-lower genes (raw and CM-adjusted). A5b transplant-only concordance of the primary program,
LVAD−TX differences and Welch p, and LVEF Spearman ρ. A6 the CS−comparator difference, Welch p and
Mann–Whitney p for all four tissue variables.

**Engine check on the A2 model.** The lme4 fits for the program score, *ACADVL* and the positive control are
singular, with patient variance 0, so `(1 | patient/section)` reduces to `(1 | section)`. The harness refits
that model with statsmodels REML, and the slopes agree. The Satterthwaite CIs (df ≈ 16,400) are read from
`a2_distance_model.tsv` and checked only for agreement with p. Because the pre-specified model has no random
slope, those degrees of freedom are spot-level for a within-section covariate. This is recorded, not
changed.

**What is not independently recomputed (taken from the producing script).** The limma fits in A1 and
A5b (reproduced to machine precision against the committed contrasts inside `a1_patient_refit.R` and
`a5b_procurement.R`). The expression-decile permutation null p values in A3 and A5b. The Hamed per-etiology
hypergeometric tests and the OLS G2−G1 fits in A5, which need the supplementary workbook
(`data/external/44161_2026_880_MOESM4_ESM.xlsx`) and are produced only by `a5_hamed.py`.

**Labelling discrepancy carried from the handoff.** The genes the handoff calls "absent annotation" are
those outside the locus-QC universe (`in_qc_universe == False`). They are not genes carrying
`flag_missing_annotation`. The ESC HF tables and Figure 1E use the "outside locus-QC universe" wording.
`flag_length_deviation` is renamed `flag_locus_quant_outlier` (the locus-quantification outlier flag) in
`results/eschf/q1_hits_annotated.tsv`.

**Specifications that could not be run as written.** A2: no CS CM-dominant spot reaches 2 mm from an
immune-dominant spot (maximum 1.73 mm), and every CS section contains immune-dominant spots. The
2–3 mm, >3 mm and no-immune-dominant-spot bins are therefore empty. A2/A6: the export contains no Kuppe normal
sections, so the normal-myocardium reference was not computed. A5: *ACTC1* is absent from the Hamed
proteome, and the sarcomeric score uses the remaining seven markers. A5b: disease duration is not in the
Neyazi supplementary tables.

### 4.1 Revision (HANDOFF_figure_fixes_ESCHF.md)

**Universe rule.** Section 0 of the handoff describes the universe as detection in CS. The committed rule in
`scripts/spatial_first/01_limma_q1q2.R` is detection fraction ≥0.25 in CS **or** in comparators, and the
handoff also asks for "the same rule as before". The committed rule was applied: 982 genes in the new
universe, 985 in the earlier one, 969 in both (`results/eschf/a0_universe.tsv`).

**Sections in the primary fit.** CS_109 has no CM-dominant spot passing the 2000-UMI floor, so the primary
fit has 11 CS sections from 8 patients (6 comparator sections from 6 patients).

**Circularity.** The program genes are hits of the primary fit, which now includes the spots nearest
immune-dominant spots. The program-score level difference in the near bins of Figure 1B is therefore not an
out-of-sample test; the OXPHOS level check (KEGG hsa00190, defined independently of the fit) is the only
non-circular level test. Whether it enters the figure, legend or text is left to Rami (handoff s.4); it is
reported in `results/eschf/a2_oxphos_level.tsv` and SUMMARY.md section 0.2.

**Mixed models.** The spot-level mixed models are kept and labelled "spot-level model (singular fits;
intervals anti-conservative)". In the primary, 3 of 4 are singular and lme4 reported a convergence warning.
The pre-specified per-patient OLS slopes (`results/eschf/a2_perpatient_slopes.tsv`) are the summary used in
Figure 1B.

**Specifications that could not be run as written.**
1. Figure 1B bins 2–3 mm, >3 mm and no-immune-dominant-spot: empty. The maximum distance from a CS
   CM-dominant spot to the nearest immune-dominant spot is 1.73 mm, and every CS section has immune-dominant
   spots. Only populated bins are drawn.
2. Universe rule: applied as committed (CS or comparators), not as worded in s.0 (CS only); see above.
3. Old-fit hard-coded count assertions in the cascade scripts (`a3a4_qc_nonfailing.py`, `a5_hamed.py`,
   `a5b_pool.py`, `a5b_procurement.R`) are skipped when `ESCHF_PRIMARY=v2`; the harness checks the new counts
   independently instead.

**Reproducibility of the earlier fit.** Every cascade script takes `ESCHF_PRIMARY` (`old` | `v2`, via
`scripts/eschf/eschf_paths.py`; R scripts read the variable directly). In `old` mode, all committed
`results/eschf/` tables are reproduced byte-for-byte (the gzip model-input file differs only in its header
timestamp).

### 4.2 Revision v3 (HANDOFF_figure_v3_ESCHF.md)

**Concordant set for Figure 1D.** The v2 concordant panel was built from the section-stable core
(`q1_core_snrna_pergene.tsv`, 18 genes), not from all spatial hits. Figure 1D now follows the handoff
(all spatial hits, same sign, snRNA FDR <5%, passing screen 3a: 24 genes) in v3; v4 reverts to the core (section 4.4). The core-based set is carried
alongside as `in_concordant_core_v2` in `results/eschf/a9_snrna_screen.tsv`. *ENG* is same sign with snRNA
FDR <5% but is not in the section-stable core, which is why it was absent from the v2 panel; it fails 3a.

**Screen 3b** is reported (`a9_screen_counts.tsv`, `a9_source_absence.tsv`) and not applied. Every 3a failure
is also flagged by 3b.

**A10 log-CPM form.** The committed patient log-CPM (`a1_patient_logcpm.tsv`) uses edgeR `cpm(prior.count = 2)`;
the plan specified the voom form. The patient pseudobulk was rebuilt exactly as `a0_primary_refit.R` builds it
(asserted to reproduce the committed log-CPM), and the voom-form log2 CPM was computed from it
(`a10_patient_log2cpm_voomform.tsv`). The KEGG fatty acid degradation and TCA sets were fetched by the same
method that reproduces the committed OXPHOS set; the A10 OXPHOS spot scores are asserted identical to the
committed ones. All three sets share genes with the primary hit list (SUMMARY.md section 5.3).

**Specifications that could not be run as written.**
1. A12: *KCTD14* is not in the snRNA count matrices; its suspect-loci row is empty.
2. A11: apex sampling cannot be recorded (Table S2 has no apex category). Subject SNU-1 is in Table S1 but has
   no snRNA sample in the GEO download (`data/external/neyazi_geo_snrna_sample_listing.txt`), and its two
   Table S2 breakdowns disagree (0 by region, 2 by histologic class).
3. A13: in the Neyazi integrated object, `individual` is a GEO sample, not a subject; samples were collapsed
   to subjects by the sample-name pattern. Several non-myeloid classes have exactly zero GPNMB detection;
   this is recorded as an observation only.
4. Foundation-model perturbation analyses were not run, as instructed.

**Table layout fix.** `a9_source_absence.tsv` now carries a `gene` index label and the `flag_3b` column for
every row (previously filled only for the FPGT reference row). No other A9 output changed on rerun.

### 4.3 A14: pathway scores with fit hits removed (HANDOFF_pathway_hits_removed.md)

Run as specified; nothing substituted. Points the handoff left open, fixed in the `a14_pathway_hits_removed.py`
docstring before running:
- Expression for decile matching is the per-gene mean voom-form log2 CPM over all 14 patients; deciles are taken
  over the 982-gene universe (not over the pool).
- The pool excludes primary-fit hits and the tested pathway's members only; genes that are hits only in the
  >=0.5 mm sensitivity fit remain in the pool, as written.
- Null percentile counts ties as half; the two-sided p uses the (1 + count) / (N + 1) form.
- The decision rule (interval excludes zero and null percentile beyond 2.5/97.5) has no direction. For OXPHOS the
  hits-removed version meets both criteria with CS higher than comparators, whereas the A10 question concerned
  lower levels in CS; SUMMARY section 6 states the criteria and the direction separately and does not interpret
  the reversal.

**Triangle list fix.** The SUMMARY Figure 1C sentence listed ten entries (each off-scale gene twice, the second
with its fill-category code printed as a value) because `a8_summary.py` selected both the `|snrna_lfc` and
`|category` rows of `a7_panel_values.tsv`. It now selects only `|snrna_lfc` rows and asserts the count equals
`n_out_of_range` (5). The figure, abstract and legend were already correct (5 triangles each, verified on the
rendered axes) and were not changed.

### 4.4 Revision v4 (user instructions, Figure 1 panels B–D)

- **Panel C.** All five screen-3a failures (*HSPB2*, *CALM2*, *NDUFC2*, *ENG*, *CFD*; these are also the only genes
  beyond the ±2.5 axis) are drawn as open neutral-grey triangles with no sign category. In v3 *NDUFC2* and *ENG* were
  filled in the "same sign, snRNA FDR <5%" style although they are excluded from the agreement fractions. The
  instruction named these two; it was applied to every screened gene, which is all five triangles.
- **Panel B legend.** The OXPHOS range sentence is replaced by the user's sentence on KEGG TCA cycle, fatty acid
  degradation and OXPHOS. The figure script asserts it against `a14_pathway_hits_removed.tsv`: for both removals
  (primary hits; primary ∪ sensitivity hits) no pathway's 95% CI lies entirely below zero. The sentence does not
  state that CS is higher for OXPHOS with primary hits removed (+0.133, 95% CI +0.019 to +0.246; section 4.3).
- **Panel D.** Uses the 18 concordant genes of the section-stable core (`in_concordant_core_v2`), replacing the
  v3 specification (all hits after 3a, 24 genes). The instruction first read "17 gene", then was corrected to "18 gene
  selection-stable core"; 18 is the section-stable core before screening. *NDUFC2* is in the core and fails 3a; it is
  kept, as the 18-gene count requires, and marked ¶ (a mark added here, not specified in the instruction). The 17-gene
  alternative is the core after 3a (`in_concordant_core_v2_screened`). Core genes not in the v3 panel: *NDUFC2* (removed by 3a in
  v3); v3 genes not in the core: *DYNLL2*, *LMO7*, *POLDIP2*, *POPDC2*, *RBPMS2*, *S100A13*, *UQCRC1*.
- Harness: `ESCHF[rev]` claims for the Panel C fill categories and the four OXPHOS legend values are replaced by
  `[v4]` claims (open-neutral drawing per triangle, Panel D set recomputed from the refits, ¶ per gene, pathway
  sentence and its six intervals).

### 4.5 Revision v5 (user instructions, Figure 1 panels C–D)

- **Panel D: 17 genes.** *NDUFC2* is removed (pooled snRNA log2FC +8.09, fails 3a, not concordant), leaving the 17
  concordant genes of the section-stable core after screen 3a (`in_concordant_core_v2_screened`). The ¶ mark and its
  key are removed.
- **snRNA vs NF column check.** Among the 18 core genes, *NDUFC2* was the only implausible value: log2FC +7.32 at mean
  expression −1.79 log2 CPM. The 17 plotted genes all lie within ±1.8 log2 (largest *IDH2* −1.79) at mean expression
  ≥ 2.62 log2 CPM (lowest *ALKBH7*); none exceeds the 3a magnitude of 3. *ALKBH7* carries a 3b flag (near-absent in
  Larson); it is not marked in D, as 3b is reported and not applied. The figure script asserts the ±3 bound. For
  context, 1,103 of 17,773 genes genome-wide exceed |log2FC| 3 in this CS vs NF contrast.
- **Panel C legend.** Added, in the user's wording: "All five screened genes were near-absent from one source study's
  quantification while present in the others, as previously observed at the FPGT-TNNI3K locus", followed by "the ±3
  log2 threshold is a pre-stated sensitivity screen, not a biological limit". The latter phrase was previously a
  parenthesis in the preceding sentence and was moved, not duplicated. Basis: all five carry the 3b flag (study mean
  < 0 log2 CPM, other four ≥ 2 higher): Reichart for *HSPB2*, *CALM2*, *NDUFC2* and *ENG*, Neyazi for *CFD*; FPGT is
  flagged in Reichart. Not captured by "present in the others": for *NDUFC2* and *CFD* a second study is also below
  0 log2 CPM (Larson, −1.75 and −1.90), although the mean of the other four is ≥ 2 higher for both.
