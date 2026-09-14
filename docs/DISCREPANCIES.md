# Where the committed analysis and the manuscript text disagree

`python scripts/manuscript_results.py` checks every quantitative claim in the brief
report against the committed tables. As of this commit, 30 claims reproduce, 8 do not,
and 2 cannot be addressed from committed data. The 8 mismatches are all in one
paragraph and all follow from one unresolved question. This document records them so
the text can be reconciled rather than the code adjusted to match it.

## 1. The published non-failing comparison (8 mismatched values)

**Manuscript text.** "Of the 139 genes lower in CS by ST, 29 were testable in DCM and
28 in HCM. Among testable genes, most were lower in disease than non-failing myocardium
(19/29 in DCM versus 47% of differential genes in the published comparison, p=0.038;
17/28 in HCM versus 56%, p=0.061)."

**What `scripts/spatial_first/06_published_reference.py` computes**, from Chaffin et al.
supplementary tables ST6 and ST7, cardiomyocyte rows, CellBender-corrected effect sizes,
background-contamination-flagged rows excluded:

| | manuscript | recomputed |
|---|---|---|
| testable in DCM | 29 | 36 |
| lower in disease, DCM | 19/29 | 25/36 |
| reference rate, DCM | 47% | 48.3% |
| binomial p, DCM | 0.038 | 0.009 |
| testable in HCM | 28 | 35 |
| lower in disease, HCM | 17/28 | 22/35 |
| reference rate, HCM | 56% | 44.7% |
| binomial p, HCM | 0.061 | 0.023 |

**The direction of the finding is unchanged and the recomputed version is stronger**:
CS-lower genes are enriched for genes lower in DCM and in HCM relative to non-failing
donors, at p = 0.009 and p = 0.023 rather than p = 0.038 and p = 0.061. The HCM
comparison moves from marginal to conventionally significant.

**What was ruled out.** The testable counts were checked against every plausible
alternative definition: gene universe (all 139 lower genes; the 115 section-stable
lower core; the QC-passing subsets; the 108 lower genes testable in both platforms;
single-flag exclusions), cell type (cardiomyocyte rows versus the pseudo-bulk rows),
effect-size column family (CellBender-corrected versus CellRanger), and the
background-contamination filter on or off. No combination yields 29 and 28 testable, and
none yields 19 and 17 lower. The discrepancy is therefore not a filter choice.

**One value does resolve.** The recomputed HCM reference rate is 44.7% lower in disease,
whose complement is 55.3% — the 56% in the text. The DCM rate in the text (47%) is
close to the recomputed 48.3%. The most likely reading is that the HCM reference rate
was taken as the fraction *not* lower in disease while the DCM rate was taken as the
fraction lower, i.e. a sign inconsistency between the two comparisons. That accounts for
the reference rates but not for the testable counts.

**Recommended resolution.** Replace the paragraph's six counts, two reference rates and
two p-values with the values in `results/spatial_first/q1_published_nf_summary.tsv`, and
state the testability caveat the script's docstring sets out: the published tables list
only genes called differential, so a CS-lower gene absent from a published table may be
unchanged in that disease or may simply not have been reported, and the binomial
reference rate is the base rate of the published differential set rather than 0.5.

## 2. Claims the committed data cannot address (2, not counted as mismatches)

**"IDH2 ... lower in all 5 genotype strata of the comparator cohort underlying our
snRNA-seq analysis."** `data/snrna/snrna_patient_cohort.tsv` carries patient, disease,
study and anatomy, not the sarcomeric/desmosomal genotype strata this claim refers to.
The stratified test was run in the single-nucleus arm against per-patient genotype
metadata that is not in this repository. Either add that metadata under `data/snrna/`
and extend `06_published_reference.py`, or drop the clause.

**"12 CS sections (11 evaluable for non-lesional myocardium)."**
`data/spatial_first/region_selection.tsv` records 12 CS sections and the post-threshold
call (7 evaluable at the 50-spot threshold), both of which reproduce. The intermediate
count of 11 — sections with any cardiomyocyte-enriched non-lesional region before the
50-spot threshold — is not a column in that table. It is consistent with the table but
not recomputed from it.

## 3. Loose wording that is numerically correct

**"Refitting the analysis by leaving out individual sections and diseases resulted in
all 192 genes retaining direction with 141 retaining significance."** All 192 retain
direction in all ten refits, so that half of the sentence covers both families. The 141
is the leave-one-*section*-out intersection (115 lower + 26 higher); the
leave-one-disease-out intersection is 117. The harness reports both. Attributing 141 to
both families overstates it slightly; splitting the sentence would be accurate.

## 4. Superseded material in the repository

The single-nucleus tables `results/DE_CS_vs_*.tsv` predate the requantification the
manuscript's cross-platform arm uses and disagree with it substantially across shared
genes. The tables the manuscript relies on are `data/snrna/snrna_CS_vs_*.tsv`. The older
directory is retained for provenance and is not on any path that produces a manuscript
number; `scripts/manuscript_results.py` does not read it.
