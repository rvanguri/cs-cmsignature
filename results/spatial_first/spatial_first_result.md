# Spatial-first analysis of normal-appearing myocardium in cardiac sarcoidosis

All thresholds in `spatial_first_config.json` were fixed before any contrast was run
(depth target grid 1500/2000/3000 UMI with 2000 as headline; detectability grid
0.10/0.25/0.40 with 0.25 as headline; 0.5 mm granuloma exclusion; 50-spot floor per section;
lineage assignment by argmax of five z-scored marker-set scores). Nothing below was
re-thresholded after seeing a result; the two places where a pre-specified choice looks
wrong in hindsight are named in section 7 and left as they were.

## 1. Cohort and region definition

21 Visium sections entered: 12 cardiac sarcoidosis (Foong, probe panel), 6 non-sarcoid
cardiomyopathies (local `new visium`: ARVC x2, Chagas, HCM, LMNA, DC-LMNA; probe panel),
4 normal donors (Kuppe, polyA). A genome-wide re-export was run on the cluster because the
prior exports carried only 342 genes.

Region definition kept spots whose dominant lineage was cardiomyocyte, after removing every
spot within 0.5 mm of an immune-dominant (myeloid or lymphoid) spot. At the headline depth
73,053 spots were scored: 30,477 cardiomyocyte-, 14,878 fibroblast-, 13,928 endothelial-,
9,789 myeloid- and 3,981 lymphoid-dominant. The argmax degeneracy flagged in the handoff is
not material: exactly 1 of 13,770 immune-dominant spots had zero counts across all five
marker sets, so immune calls are driven by signal rather than by tie-breaks on empty spots.

**The exclusion is severe and asymmetric.** 5 of 12 CS sections fall below the 50-spot
floor (CS_103 and CS_109 retain 0 spots, CS_104 retains 6, CS_106 22, CS_101-1 30), leaving
7 CS, 6 comparator and 4 normal sections evaluable at the headline depth (Fig 1a,b). All 6
comparators and all 4 normals clear the floor. Because power is set by the smaller arm, the
whole analysis rests on 7 CS sections.

A caveat on the proxy itself: the immune-dominant fraction is *highest in the polyA normals*
(25.9% of spots, vs 22.5% in CS and 9.0% in comparators). The argmax proxy therefore tracks
assay chemistry as well as lesion content, and the "granuloma exclusion" is better described
as exclusion around immune-marker-high spots.

## 2. Gene space

The three sources do not share a reference. The probe panels carry 18,085 (CS) and 18,132
(comparator) symbols; the polyA normals carry 35,468. The inner join used for every contrast
is **17,619 symbols** — the normals lose 17,849 symbols to the join, the CS sections 466 and
the comparators 513 (Fig 2a). Per-section accounting is in `gene_space_accounting.tsv`.

## 3. Q2 — disease vs normal donor: not interpretable as disease biology

At the headline cell (T = 2000, detectability 0.25) 1,026 genes were testable and **694 of
them (68%) are differential at FDR 5%**, median |log2FC| 1.20. Three independent features
identify this as an assay-chemistry contrast:

- **Mitochondrial read fraction splits by source, not by disease**: median 38.4% in the polyA
  normals vs 11.5% (CS) and 7.8% (comparators), a 3-to-5-fold difference (Fig 2b). Dropping
  all MT genes barely changes the outcome (659 of 1,015 genes still differential), so the MT
  genes are a symptom of the chemistry difference rather than its only carrier.
- **Three genes are detected in normals only**, hard presence/absence calls rather than fold
  changes (open symbols, Fig 2c); they carry `flag_source_specific_absence` in
  `q2_locus_qc.tsv`.
- **The two probe-panel cohorts shift identically against the normals**: CS-vs-normal and
  comparator-vs-normal effects correlate at r = 0.83 along the identity line (Fig 2d). A
  disease effect would not be shared this closely by six unrelated cardiomyopathies.

The pre-specified Q2 disease-enriched list is delivered as `q2_raw_T2000_d25.tsv` plus
`q2_locus_qc.tsv` **with the confound attached, and no target shortlist is nominated from
it.** 305 of the 694 hits replicate in direction across both disease sources, but
replication across two cohorts sharing one chemistry against normals using another does not
address the confound. Answering Q2 requires polyA CS sections or probe-panel normals;
neither exists in this data.

## 4. Q1 — CS vs other cardiomyopathy, within one chemistry: not null

Both arms of Q1 are probe-panel data, so the section 3 confound does not apply. Of 985
testable genes, **192 are differential at FDR 5% (53 up in CS,
139 down), median |log2FC| 0.91**; 128 pass locus QC
(Fig 3a). Flags among the 192: 31 unannotated, 20 in segmental duplications, 14 with a
gene-span deviation, 0 presence/absence, 0 readthrough.

Three controls, all in `q1_control_comparison.tsv` (176 of the 192 hits testable in all
three):

| control | question | result |
|---|---|---|
| Panel v2 vs v1 among comparators | is the signal a probe-panel version effect? | r = +0.01 with the Q1 axis; 3/176 hits differential in the control; 25 of 895 genes genome-wide |
| Panel-matched comparators only | does it survive dropping chemistry-mismatched comparators? | r = +0.97, 176/176 keep sign, 74/176 still FDR<5% (117 of 985 genome-wide) |
| Adjusted for cardiomyocyte content | is it a purity shift created by the CS-only lesion exclusion? | r = +0.99, 138/176 still FDR<5% (169 of 985 genome-wide) |

**The bounded-null claim in the handoff is refuted for effects of this size.** The bound
itself: with 7 CS vs 6 comparator sections, median posterior residual SD 0.79
and 15.1 residual df, the standard error on a log2 fold change is 0.44, giving 80%
power at |log2FC| = 1.32 against a nominal 0.05 threshold and 2.85 against a
multiplicity-corrected one. Effects smaller than roughly one doubling are not addressable in
this cohort, and the 192 hits sit at median 0.91 — they are detected through the moderated
variance model, not by exceeding the naive per-gene bound.

## 5. What the CS-enriched genes are

The 40 up-in-CS genes passing locus QC are dominated by respiratory-chain and
ATP-synthase subunits (NDUFS7, NDUFB7, COX17, ATP5MC3, ATP5MG), sarcomeric and cytoskeletal
genes (ANKRD2, MYO18B, NRAP, MYOM2, TNNI3, LDB3), and cytosolic metabolic enzymes (MDH1,
LPL, MPC1, ECI1). Full list: ANKRD2, EEF2, MYO18B, ATP5MC3, NRAP, CAVIN2, MDH1, MYOM2, COX17, ADPRHL1, LPL, NDUFS7, ALPK3, PRNP, TNNI3, LAMC1, MYOM3, LGALS1, CPVL, SRL, NDUFB7, MAP7D1, ATP5MG, LDB3, PTGFRN, PDLIM1, TMX4, TSPYL1, TP53INP2, ECI1, COX6A2, MPC1, CLIC5, PSMC5, LMO7, DYNLL2, PARM1, RNF10, TRDN, ATP2B4.

This is a cardiomyocyte oxidative and structural program, **not an immune program** — no
cytokine, chemokine, complement or antigen-presentation gene appears in it. Read literally,
the normal-appearing myocardium of CS hearts looks metabolically and sarcomerically *more
intact* than the myocardium of the six comparator cardiomyopathies. Two readings survive the
controls in section 4 and this analysis cannot separate them: (i) genuine relative
preservation of the myocyte program in non-lesional CS myocardium, or (ii) residual
compositional structure that a scalar myocyte-fraction covariate does not fully absorb,
since the lesion-adjacency exclusion hit CS sections far harder than comparators. The
purity-adjusted contrast argues against (ii) being the whole story but does not exclude it.

Target mapping (Open Targets, 40 genes queried by Ensembl ID) is reported as
database provenance, not as a nomination: 7 genes carry approved drugs or clinical
candidates, 13 carry small-molecule tractability evidence, 0 carry chemical probes. Given
the character of the list, nominating any of these as a cardiac-sarcoidosis target would be
over-reading a myocyte-metabolism axis.

| gene | log2FC | FDR | target class | drugs/candidates | small-molecule tractability |
|---|---|---|---|---|---|
| ANKRD2 | +2.00 | 0.022 | — | 0 | — |
| EEF2 | +0.94 | 0.008 | Other cytosolic protein | 10 | Structure with Ligand |
| MYO18B | +0.79 | 0.010 | — | 0 | — |
| ATP5MC3 | +0.77 | 0.031 | — | 0 | — |
| NRAP | +0.76 | 0.002 | — | 0 | — |
| CAVIN2 | +0.75 | 0.048 | — | 0 | — |
| MDH1 | +0.75 | 0.029 | Enzyme;Oxidoreductase | 0 | — |
| MYOM2 | +0.71 | 0.020 | — | 0 | — |
| COX17 | +0.70 | 0.012 | — | 0 | — |
| ADPRHL1 | +0.68 | 0.009 | — | 0 | — |
| LPL | +0.68 | 0.012 | Enzyme;Hydrolase | 0 | Structure with Ligand;High-Quality Ligand;Dr |
| NDUFS7 | +0.66 | 0.040 | Enzyme;Oxidoreductase | 3 | Approved Drug;Structure with Ligand |
| ALPK3 | +0.66 | 0.006 | Enzyme;Transferase | 0 | — |
| PRNP | +0.64 | 0.009 | Surface antigen | 0 | Med-Quality Pocket |
| TNNI3 | +0.63 | 0.024 | — | 1 | Advanced Clinical;Structure with Ligand;High |

Full table: `q1_target_mapping_opentargets.tsv`.

## 6. Sensitivity grid and the confirmatory single-nucleus arm

Both arms behave the same way in all 9 grid cells: Q2 declares 54-77% of tested genes
differential, Q1 declares 18-23%, with stable median effect sizes. Nothing about the
headline cell is special.

**Q2 (disease vs normal)**

| depth | detect. | tested | differential | median abs log2FC | evaluable sections |
|---|---|---|---|---|---|
| 1500 | 0.10 | 2503 | 1348 (54%) | 1.17 | 16 |
| 1500 | 0.25 | 720 | 501 (70%) | 1.27 | 16 |
| 1500 | 0.40 | 363 | 280 (77%) | 1.36 | 16 |
| 2000 | 0.10 | 3337 | 1845 (55%) | 1.1 | 17 |
| 2000 | 0.25 | 1026 | 694 (68%) | 1.2 | 17 |
| 2000 | 0.40 | 510 | 395 (77%) | 1.26 | 17 |
| 3000 | 0.10 | 4692 | 2555 (54%) | 1.04 | 17 |
| 3000 | 0.25 | 1688 | 1084 (64%) | 1.1 | 17 |
| 3000 | 0.40 | 817 | 583 (71%) | 1.18 | 17 |

**Q1 (CS vs other cardiomyopathy)**

| depth | detect. | tested | differential | up in CS | median abs log2FC |
|---|---|---|---|---|---|
| 1500 | 0.10 | 2483 | 483 (19%) | 151 | 0.94 |
| 1500 | 0.25 | 694 | 130 (19%) | 37 | 0.88 |
| 1500 | 0.40 | 331 | 60 (18%) | 14 | 0.9 |
| 2000 | 0.10 | 3291 | 756 (23%) | 229 | 0.92 |
| 2000 | 0.25 | 985 | 192 (19%) | 53 | 0.91 |
| 2000 | 0.40 | 479 | 89 (19%) | 21 | 0.98 |
| 3000 | 0.10 | 4689 | 1099 (23%) | 326 | 0.88 |
| 3000 | 0.25 | 1681 | 384 (23%) | 116 | 0.89 |
| 3000 | 0.40 | 788 | 157 (20%) | 39 | 0.9 |

Evaluable section counts are 17 at T = 2000 and 3000 but 16 at T = 1500, where one CS section
drops out — the opposite of the naive expectation. A lower depth target retains more sparse
spots, some of which are called immune-dominant, which enlarges the 0.5 mm excluded
neighbourhood and shrinks the cardiomyocyte region.

The single-nucleus arm cannot confirm Q1, for two separate reasons (Fig 4):

- **The CS-vs-HCM contrast is a gene-length artifact.** Effect size regresses on gene span
  with slope +1.43 log2FC per decade (r = +0.34, P below double precision), and 10,580 of
  14,446 genes are called differential. No other contrast behaves this way: CS-vs-DCM
  -0.15, CS-vs-ARVC -0.33, CS-vs-NF +0.06. This closes the handoff's standing suspicion
  about that contrast as an artifact.
- **The two assays share no effect axis.** Across 922 genes testable in both, the spatial
  Q1 axis correlates with the single-nucleus CS-vs-other axes at r = -0.05 (HCM), +0.02
  (DCM), +0.01 (ARVC); of 49 spatial up-in-CS genes, 33-59% agree in direction, i.e. chance.
  The single-nucleus data neither supports nor contradicts the spatial Q1 result.

## 7. Thresholds that look wrong in hindsight (reported, not changed)

1. **0.5 mm granuloma exclusion** costs 5 of 12 CS sections, and because the immune-marker
   proxy is itself chemistry-dependent (section 1) it removes tissue asymmetrically across
   cohorts. A smaller radius, or a lesion mask drawn on histology rather than on an
   expression argmax, would keep the CS arm at 10-12 sections.
2. **Detectability 0.25** shrinks the testable space from 3,291 genes (at 0.10) to 985, a
   3.3-fold loss of gene space to guard against sparse-count instability that the voom
   weights already model.
3. The pre-specified design has no route to a chemistry-free Q2 contrast. That was a design
   defect that should have been caught before the contrast was pre-specified.

## 8. Bottom line

- Q2 as specified is an assay-chemistry contrast and yields no disease-enriched target list.
- Q1 is not null: 192 genes at FDR 5% within one chemistry, robust to probe-panel version,
  to dropping chemistry-mismatched comparators, and to myocyte-content adjustment.
- The CS-enriched half of that signature is a cardiomyocyte oxidative and sarcomeric
  program, not an immune one, and it is not a target list.
- The cohort bounds detectable effects at roughly |log2FC| 1.3; anything subtler is untested
  here, on 7 CS sections after exclusion.

## 9. Follow-up: the dropped CS sections, and Q1 without mitochondrial genes

### 9.1 Scoring the 192-gene axis in the sections the exclusion dropped

The axis is the signed mean z-score of the 192 Q1 hits (weights = sign of log2FC; 53 up,
139 down), with per-gene mean and SD taken from the 13 evaluable probe-panel sections only,
so the five dropped CS sections are scored out-of-sample. Scores are in
`q1_axis_scores.tsv`, the depth control in `q1_axis_depth_control.tsv` (Fig 5).

**Two of the five cannot be scored at all**: CS_103 and CS_109 retain zero cardiomyocyte
spots after the 0.5 mm exclusion, so they have no pseudobulk under the pre-specified region.
Scoring them would require re-exporting those sections under a different region definition.

The other three score on the CS side of the axis, out-of-sample:

| section | spots | UMI | axis score | permutation z | axis genes detected |
|---|---|---|---|---|---|
| CS_104 | 6 | 12k | +2.09 | +2.6 | 131/192 |
| CS_106 | 21 | 42k | +1.17 | +5.3 | 171/192 |
| CS_101-1 | 29 | 58k | +0.99 | +8.6 | 177/192 |

For reference, the 7 evaluable CS sections span -0.51 to +1.13 and the 6 comparators span
-0.73 to -0.63; separation is complete on the training sections, though that is in-sample
and therefore circular. Two evaluable CS sections (CS_102-2, CS_101-2) sit at negative
absolute scores while still ranking above every comparator.

Because 139 of the 192 genes are down-in-CS, dropout at low depth pushes the score upward
mechanically, so the raw values above are not directly comparable to the full-depth
sections. Multinomial downsampling of the 13 training sections to 12k / 42k / 58k UMI
(20 replicates each) bounds that bias: comparator sections stay negative at every matched
depth (mean -0.57 at 12k, drift of only
+0.11 from full depth), while evaluable CS
sections inflate by +0.95 at 12k and by +0.14 to +0.21 at 42-58k. So CS_104's +2.09 is
substantially depth-inflated and its depth-matched CS reference range is wide
(-0.73 to +2.87), whereas CS_106 and CS_101-1 sit inside the depth-matched CS range and far
outside the comparator range.

**Interpretation.** The exclusion did not select a CS subgroup that happens to carry the
signature: the sections it dropped — the most lesion-burdened ones — carry the same axis,
scored out-of-sample and against a depth-matched control. This addresses a selection
concern but not the purity concern, since the dropped sections are also the ones whose
retained tissue is smallest and most lesion-adjacent.

### 9.2 Q1 with mitochondrial genes removed and counts renormalized

All 11 MT genes were dropped from the pseudobulk matrix and TMM factors and voom weights
recomputed from the MT-free counts (the MT genes carry 10.8% of the CS library and 7.8% of
the comparator library, so this is a small renormalization within the probe-panel arm).

| run | tested | FDR<5% | up in CS | median abs log2FC |
|---|---|---|---|---|
| original | 985 | 192 | 53 | 0.908 |
| MT-free, renormalized | 974 | 197 | 60 | 0.871 |

Across the 974 jointly tested genes the effect estimates are unchanged to four decimal
places in correlation (r = +1.0000); 190 of the 192 original hits remain, 2 drop out
(AGTPBP1, BCL2L13, both at |log2FC| ~ 0.5) and 7 enter (CUX1, CYB5R1, ESRRA, MYOM1, PRDX5,
SGCA, UBXN6, median |log2FC| 0.46). Q1 is not a mitochondrial-content effect — unlike Q2,
where MT genes are one visible symptom of a chemistry difference that persists after their
removal. Per-gene comparison: `q1_noMT_comparison.tsv`.

## 10. Follow-up: leave-one-comparator-section-out jackknife on Q1

Six comparator sections cover five diseases (ARVC x2, Chagas, LMNA_DCM x2, HCM), so a
section-specific contrast could masquerade as a disease one. Q1 was refit six times, each
time marking one comparator section non-evaluable, with the gene space held fixed at the
same 985 detectability-passing genes so runs are directly comparable
(`q1_leave_one_comparator_out.tsv`, Fig 6).

| dropped | disease | FDR<5% | up in CS | r (985 genes) | r (192 hits) | sign concordant | retained of 192 | new hits |
|---|---|---|---|---|---|---|---|---|
| − LMNA | LMNA_DCM | 151 | 31 | 0.9971 | 0.9988 | 192/192 | 148 | 3 |
| − DC_LMNA | LMNA_DCM | 166 | 40 | 0.9970 | 0.9981 | 192/192 | 157 | 9 |
| − ARVC_2 | ARVC | 178 | 49 | 0.9972 | 0.9987 | 192/192 | 169 | 9 |
| − HCM | HCM | 186 | 50 | 0.9913 | 0.9987 | 192/192 | 167 | 19 |
| − ARVC_1 | ARVC | 191 | 51 | 0.9949 | 0.9978 | 192/192 | 174 | 17 |
| − Chagas | Chagas | 192 | 53 | 0.9946 | 0.9984 | 192/192 | 172 | 20 |

The result does not collapse on any single section. Hit counts range 151-192 against the
full-cohort 192; the effect vector correlates with the full-cohort axis at r >= 0.991 over
all 985 tested genes and r >= 0.998 over the 192 hits; and sign concordance on the 192 is
192/192 in every one of the six runs. Genes that fall out do so by crossing the FDR line
from just below it, not by reversing (Fig 6c): in the worst run the 44 lost genes all sit
within a factor of ~2 of the threshold. **141 of the 192 are significant in the full cohort
and in all six leave-one-out runs** (26 of them up in CS) — that intersection is the
defensible core of the axis.

The sensitive quantity is not the axis but the up-in-CS count, which falls from 53 to 31
when either LMNA_DCM section is dropped. The up direction is already the one carrying the
cardiomyocyte-purity concern from section 6, and this jackknife adds a second reason to
treat the up-in-CS list as the weaker half of the result. The down direction is stable
across all six runs.

Caveat on interpretation: with two sections per disease for ARVC and LMNA_DCM, dropping one
of a pair leaves that disease still represented, so this jackknife bounds section-level
influence, not disease-level influence. A leave-one-disease-out refit would drop to four
diseases and is a different, harsher test that was not run.

## 11. Follow-up: core-set concordance against a pooled snRNA comparator contrast

The two jackknife-stable halves of the Q1 axis were tested separately against the
single-nucleus data. The snRNA comparator contrast is CS vs DCM and CS vs ARVC combined by
fixed-effect inverse-variance meta-analysis (per-gene SE recovered as |log2FC / t|);
**CS vs HCM is excluded** because that contrast carries a strong gene-length bias
(slope +1.43 log2FC per log10 span, section 8) and CS vs NF is a disease-vs-normal contrast,
not a comparator one. 17,747 genes pool; the null used throughout is 2000 draws of
expression-matched genes (AveExpr deciles), which absorbs any global sign asymmetry in the
snRNA contrast.

| set | in snRNA | sign-concordant | matched null | perm p | replicating FDR<5% | null rep. rate | perm p | Spearman rho | p |
|---|---|---|---|---|---|---|---|---|---|
| stable down (115) | 108 | 65 (60.2%) | 47.6% | 0.006 | 15 (13.9%) | 12.1% | 0.32 | +0.06 | 0.53 |
| stable up (26) | 26 | 12 (46.2%) | 52.0% | 0.79 | 3 (11.5%) | 13.0% | 0.69 | **-0.51** | 0.008 |

**The two halves behave differently.** The 115 stable down genes carry a real but weak sign
lean: 60% concordant against a 48% expression-matched null (p = 0.006), which is a
directional signal and nothing more — gene-level replication (13.9% at pooled FDR 5%) is
indistinguishable from the matched background rate (12.1%, p = 0.32), and effect magnitudes
do not correspond at all (Spearman rho = +0.06). The 15 that do replicate sign-matched are
ACADVL, IDH2, ALKBH7, AURKAIP1, POPDC2, SLC2A4, FBXW5, JUP, TINAGL1, COL1A1, SPTAN1, FITM2,
SPAG7, ARL6IP4, MPRIP.

The 26 stable up genes show no concordance at all (46% vs a 52% null, p = 0.79) and,
more informatively, an **inverted** magnitude relationship: Spearman rho = -0.51 (p = 0.008),
meaning the genes with the largest spatial up-in-CS effect are the ones most negative in the
pooled snRNA contrast. The four largest spatial effects — ANKRD2 (+2.00), EEF2 (+0.94),
CITED2 (+0.80), MYO18B (+0.79) — all go the other way in snRNA, EEF2 significantly
(pooled log2FC -0.86, FDR 0.001). Only DYNLL2, LMO7 and TRDN replicate sign-matched, which
is at the background rate for a 26-gene set.

This is the third independent line pointing the same way: the up-in-CS half of Q1 behaves
like a cardiomyocyte-content artifact of the lesion exclusion (section 6 purity adjustment,
section 10 jackknife instability, and now sign-inversion in nucleus-resolved data where
per-cell identity is assigned rather than inferred from spot mixture), while the down-in-CS
half is stable and weakly reproducible in direction. **No target should be nominated from
the up-in-CS list.**

Caveat: the two pooled contrasts share the same CS samples, so the pooled SE understates
uncertainty and `fdr_pool` is anti-conservative. The sign-concordance and permutation
results above do not depend on the pooled SE; only the FDR<5% replication columns do, and
those are the ones reported as not exceeding background. Per-gene values:
`q1_core_snrna_pergene.tsv`; summary: `q1_core_snrna_concordance.tsv` (Fig 7).

## 12. Follow-up: leave-one-disease-out jackknife

The comparator arm covers four diseases, two of them contributed by two sections each, so the
section-level jackknife of section 10 bounds section influence but not disease influence.
Q1 was refit dropping each disease entirely — both ARVC sections together, both LMNA_DCM
sections together, and the single Chagas and HCM sections (those two runs are identical to
their section-level counterparts and are reused, not refit). Gene space is again held fixed
at the same 985 genes (`q1_leave_one_disease_out.tsv`).

| dropped disease | sections dropped | comparator sections left | FDR<5% | up in CS | r (985 genes) | r (192 hits) | sign concordant | retained of 192 | new hits |
|---|---|---|---|---|---|---|---|---|---|
| − LMNA_DCM | 2 | 4 | 134 | 29 | 0.9886 | 0.9955 | 192/192 | 124 | 10 |
| − ARVC | 2 | 4 | 170 | 45 | 0.9867 | 0.9949 | 192/192 | 148 | 22 |
| − HCM | 1 | 5 | 186 | 50 | 0.9913 | 0.9987 | 192/192 | 167 | 19 |
| − Chagas | 1 | 5 | 192 | 53 | 0.9946 | 0.9984 | 192/192 | 172 | 20 |

The harsher test gives the same verdict. Dropping a whole disease costs more power — n falls
from 13 to 11 sections when a paired disease goes — and the hit count falls further, to 134
when LMNA_DCM is removed entirely and 170 when ARVC is. But the axis itself does not move:
r >= 0.987 over all 985 tested genes, r >= 0.995 over the 192 hits, and sign concordance is
again 192/192 in every run. LMNA_DCM is the most influential disease and ARVC second, the
same ordering as at section level, and again it is the up-in-CS count that absorbs most of
the loss (53 to 29 without LMNA_DCM).

**A disease-level stable core of 115 genes** (95 down, 20 up in CS) is significant in the
full cohort and in all four disease-drop runs. It is a strict subset of the 141-gene
section-level core: 95 of the 115 section-stable down genes and 20 of the 26 section-stable
up genes survive. Note that this 115 is a different set from the 115 section-stable down
genes listed below — the coincidence of size is accidental.

Caveat: with four diseases, each run leaves three, and no jackknife on four units can
distinguish "robust across diseases" from "robust across the three that remain". This
bounds single-disease influence and nothing stronger.

## 13. The 115 section-stable down-in-CS genes

Genes down in CS relative to the comparator arm at FDR 5% in the full cohort and in all six
section-level leave-one-out runs (section 10), ordered by effect size. `worst section-LOO
FDR` is the largest FDR the gene took across those six runs; `survives disease-LOO` marks
the 95 that are also significant in all four disease-level runs (section 12); the snRNA
columns are the pooled DCM + ARVC comparator contrast of section 11, where concordant means
the pooled effect is also negative. Machine-readable: `q1_stable_down_115.tsv`.

| gene | log2FC | FDR | worst section-LOO FDR | survives disease-LOO | pooled snRNA log2FC | snRNA FDR | concordant |
|---|---|---|---|---|---|---|---|
| ATP5F1A | -6.09 | 0.015 | 0.029 | yes | — | — | - |
| SPTAN1 | -5.80 | 0.005 | 0.011 | yes | -0.25 | 0.014 | yes |
| TJP1 | -5.61 | 0.006 | 0.012 | yes | -0.05 | 0.743 | yes |
| PLEC | -5.25 | 0.006 | 0.012 | yes | -0.25 | 0.100 | yes |
| TIMP3 | -5.15 | 0.000842 | 0.000949 | yes | -0.02 | 0.973 | yes |
| IDH2 | -4.42 | 0.000359 | 0.000455 | yes | -0.59 | 0.008 | yes |
| TIMM8B | -4.35 | 0.018 | 0.020 | yes | -0.04 | 0.908 | yes |
| SAMD4A | -4.26 | 0.024 | 0.026 | yes | -0.18 | 0.226 | yes |
| ASPH | -4.12 | 0.001 | 0.000955 | yes | 0.15 | 0.038 | no |
| SPAG7 | -4.09 | 0.017 | 0.016 | yes | -0.47 | 0.025 | yes |
| SOD1 | -3.85 | 0.023 | 0.022 | yes | -0.32 | 0.238 | yes |
| TINAGL1 | -3.84 | 0.008 | 0.006 | yes | -0.92 | 1.03e-11 | yes |
| TMEM59 | -3.83 | 0.024 | 0.041 | yes | -0.23 | 0.283 | yes |
| SLC4A3 | -3.69 | 5.05e-06 | 1.66e-05 | yes | -0.11 | 0.390 | yes |
| CSDE1 | -3.57 | 0.003 | 0.003 | yes | 0.16 | 0.293 | no |
| TMEM14C | -3.31 | 0.009 | 0.009 | yes | -0.18 | 0.435 | yes |
| ITGA7 | -3.16 | 8.26e-05 | 0.000116 | yes | -0.01 | 0.934 | yes |
| PIK3R1 | -3.08 | 0.000123 | 0.000234 | yes | 0.24 | 0.317 | no |
| FKBP5 | -3.05 | 0.012 | 0.013 | yes | 0.37 | 0.436 | no |
| ATRX | -3.04 | 0.003 | 0.003 | yes | 0.08 | 0.215 | no |
| TM9SF2 | -3.02 | 0.022 | 0.020 | yes | 0.09 | 0.240 | no |
| THBS4 | -3.00 | 0.003 | 0.003 | yes | 0.12 | 0.751 | no |
| FH | -2.92 | 0.002 | 0.002 | yes | -0.04 | 0.834 | yes |
| TIMMDC1 | -2.81 | 7.67e-05 | 0.00011 | yes | 0.07 | 0.701 | no |
| HNRNPH3 | -2.46 | 8.62e-05 | 0.000175 | yes | 0.24 | 0.002 | no |
| ITGAV | -2.44 | 4.13e-05 | 0.00011 | yes | 0.01 | 0.958 | no |
| NEXN | -2.27 | 0.008 | 0.017 | yes | 0.12 | 0.396 | no |
| HSPB7 | -2.22 | 0.010 | 0.023 | yes | -0.05 | 0.799 | yes |
| KLHDC2 | -2.17 | 4.13e-05 | 0.00011 | yes | 0.01 | 0.972 | no |
| HRC | -2.16 | 0.003 | 0.007 | yes | 0.29 | 0.285 | no |
| JUP | -2.14 | 0.027 | 0.025 | yes | -0.42 | 0.000134 | yes |
| HSD17B4 | -1.86 | 0.000498 | 0.000635 | yes | -0.13 | 0.233 | yes |
| FBXO32 | -1.85 | 0.001 | 0.002 | yes | -0.33 | 0.181 | yes |
| HSPB8 | -1.83 | 0.006 | 0.012 | yes | 0.03 | 0.872 | no |
| FN1 | -1.77 | 0.003 | 0.007 | yes | -2.01 | 0.141 | yes |
| HECTD1 | -1.77 | 0.000621 | 0.00073 | yes | -0.16 | 0.078 | yes |
| COX6B1 | -1.76 | 0.003 | 0.005 | yes | -0.45 | 0.124 | yes |
| BAG6 | -1.76 | 4.13e-05 | 0.00011 | yes | -0.21 | 0.085 | yes |
| KIDINS220 | -1.76 | 4.13e-05 | 0.00011 | yes | 0.16 | 0.045 | no |
| ANK3 | -1.74 | 0.007 | 0.004 | yes | -0.03 | 0.888 | yes |
| ARL6IP4 | -1.71 | 0.000123 | 0.000506 | yes | -0.38 | 0.001 | yes |
| STAU2 | -1.70 | 0.000123 | 0.000586 | yes | 0.96 | 1.72e-35 | no |
| AK1 | -1.70 | 0.000842 | 0.003 | yes | 1.01 | 2.74e-05 | no |
| PICALM | -1.69 | 0.001 | 0.002 | yes | -0.07 | 0.724 | yes |
| QKI | -1.65 | 0.003 | 0.006 | yes | 0.15 | 0.283 | no |
| AURKAIP1 | -1.62 | 0.002 | 0.007 | yes | -0.62 | 0.004 | yes |
| ATP5MC1 | -1.62 | 0.006 | 0.005 | yes | — | — | - |
| IMMT | -1.59 | 8.62e-05 | 0.000691 | yes | -0.00 | 0.993 | yes |
| CAPZA2 | -1.57 | 0.000442 | 0.00073 | yes | -0.02 | 0.848 | yes |
| KPNA4 | -1.55 | 0.038 | 0.041 | yes | 0.17 | 0.360 | no |
| FNIP2 | -1.54 | 0.005 | 0.010 | yes | 0.06 | 0.838 | no |
| STK38L | -1.46 | 0.000581 | 0.002 | yes | 0.16 | 0.389 | no |
| KANK1 | -1.44 | 4.13e-05 | 0.00011 | yes | -0.02 | 0.913 | yes |
| COL1A1 | -1.42 | 0.024 | 0.038 | yes | -2.61 | 0.005 | yes |
| FNDC3B | -1.40 | 0.002 | 0.005 | yes | 0.05 | 0.800 | no |
| PLAAT3 | -1.40 | 0.002 | 0.003 | yes | — | — | - |
| FGF12 | -1.36 | 0.009 | 0.026 | no | -0.14 | 0.716 | yes |
| ATP5PD | -1.36 | 0.013 | 0.026 | yes | — | — | - |
| HMGN3 | -1.31 | 0.001 | 0.002 | yes | -0.25 | 0.218 | yes |
| KIF1B | -1.30 | 0.003 | 0.009 | yes | -0.15 | 0.185 | yes |
| HSF1 | -1.30 | 0.001 | 0.003 | yes | -0.14 | 0.346 | yes |
| FITM2 | -1.29 | 0.029 | 0.038 | yes | -0.38 | 5.3e-05 | yes |
| BAG1 | -1.28 | 0.013 | 0.029 | no | -0.02 | 0.906 | yes |
| STRAP | -1.27 | 0.000442 | 0.001 | yes | -0.24 | 0.108 | yes |
| GNB1 | -1.26 | 0.009 | 0.022 | yes | -0.07 | 0.772 | yes |
| FOXO3 | -1.26 | 0.006 | 0.014 | yes | 0.01 | 0.972 | no |
| HSPB2 | -1.24 | 0.000258 | 0.001 | yes | 7.05 | 1.13e-58 | no |
| ACTN1 | -1.24 | 0.002 | 0.004 | yes | -0.20 | 0.358 | yes |
| ALKBH7 | -1.20 | 0.001 | 0.002 | yes | -0.92 | 0.00014 | yes |
| STAU1 | -1.18 | 0.000359 | 0.001 | yes | -0.02 | 0.880 | yes |
| RNF115 | -1.17 | 0.015 | 0.034 | no | -0.07 | 0.669 | yes |
| USP9X | -1.12 | 0.000797 | 0.002 | yes | -0.02 | 0.880 | yes |
| LRRC39 | -1.10 | 0.000842 | 0.002 | yes | -0.16 | 0.432 | yes |
| AK4 | -1.06 | 0.002 | 0.003 | yes | 0.35 | 0.013 | no |
| FBXW5 | -1.03 | 0.015 | 0.025 | yes | -0.34 | 0.020 | yes |
| ATP5MC2 | -1.00 | 0.017 | 0.043 | no | — | — | - |
| KANK2 | -0.97 | 0.002 | 0.006 | yes | -0.04 | 0.801 | yes |
| IGFBP2 | -0.96 | 0.009 | 0.020 | yes | 0.13 | 0.772 | no |
| SYNPO | -0.96 | 0.012 | 0.028 | no | -0.15 | 0.516 | yes |
| SPTB | -0.94 | 0.013 | 0.028 | no | 0.24 | 0.019 | no |
| ACADVL | -0.94 | 0.006 | 0.012 | yes | -0.67 | 0.000368 | yes |
| COX7A2 | -0.92 | 0.008 | 0.016 | yes | 0.11 | 0.724 | no |
| NRDC | -0.91 | 0.001 | 0.004 | yes | -0.06 | 0.474 | yes |
| NDUFAF8 | -0.91 | 0.006 | 0.012 | yes | — | — | - |
| ARL6IP5 | -0.90 | 0.006 | 0.017 | yes | -0.22 | 0.358 | yes |
| ACADM | -0.89 | 0.022 | 0.048 | no | 0.10 | 0.686 | no |
| USP53 | -0.88 | 0.005 | 0.011 | yes | 0.14 | 0.280 | no |
| MEF2A | -0.84 | 0.002 | 0.006 | yes | 0.29 | 0.004 | no |
| SOD2 | -0.83 | 0.003 | 0.006 | yes | — | — | - |
| IDH3B | -0.83 | 0.008 | 0.022 | yes | 0.74 | 1.2e-05 | no |
| HK1 | -0.82 | 0.002 | 0.005 | yes | 0.03 | 0.871 | no |
| RBM33 | -0.75 | 0.010 | 0.027 | no | 0.16 | 0.051 | no |
| NCKAP1 | -0.75 | 0.003 | 0.011 | yes | 0.09 | 0.250 | no |
| MPRIP | -0.72 | 0.001 | 0.004 | yes | -0.43 | 0.000676 | yes |
| NEK7 | -0.71 | 0.024 | 0.047 | no | 0.09 | 0.565 | no |
| PAM | -0.70 | 0.017 | 0.034 | no | -0.10 | 0.655 | yes |
| NCOA4 | -0.68 | 0.013 | 0.026 | no | 0.25 | 0.168 | no |
| TOMM7 | -0.67 | 0.022 | 0.045 | no | -0.48 | 0.102 | yes |
| SDHB | -0.64 | 0.006 | 0.014 | yes | -0.22 | 0.212 | yes |
| ANK2 | -0.63 | 0.008 | 0.020 | yes | -0.08 | 0.796 | yes |
| PFKP | -0.61 | 0.022 | 0.044 | no | 0.17 | 0.258 | no |
| CAV2 | -0.59 | 0.018 | 0.033 | no | 0.52 | 9.05e-05 | no |
| PNRC1 | -0.57 | 0.017 | 0.037 | no | 0.37 | 0.046 | no |
| NPTN | -0.57 | 0.010 | 0.025 | yes | -0.09 | 0.460 | yes |
| EIF5 | -0.56 | 0.018 | 0.034 | yes | -0.20 | 0.146 | yes |
| PTPRM | -0.56 | 0.013 | 0.028 | no | 0.08 | 0.783 | no |
| HDLBP | -0.56 | 0.017 | 0.028 | yes | 0.06 | 0.770 | no |
| AIFM1 | -0.55 | 0.013 | 0.025 | yes | -0.12 | 0.437 | yes |
| ATPAF1 | -0.54 | 0.018 | 0.038 | no | -0.01 | 0.978 | yes |
| EMC4 | -0.53 | 0.022 | 0.035 | yes | -0.12 | 0.394 | yes |
| SLC2A4 | -0.52 | 0.024 | 0.044 | yes | -0.96 | 2.26e-09 | yes |
| POPDC2 | -0.48 | 0.022 | 0.041 | no | -0.28 | 0.020 | yes |
| SYNPO2 | -0.46 | 0.010 | 0.025 | yes | 0.15 | 0.401 | no |
| PDHB | -0.43 | 0.024 | 0.048 | no | -0.33 | 0.074 | yes |
| VEGFA | -0.38 | 0.024 | 0.044 | no | -0.16 | 0.341 | yes |
