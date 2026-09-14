# Verification: what reproduces, and what it cannot

`python scripts/manuscript_results.py` checks every quantitative claim in the brief
report against the committed tables. As of this commit, **40 claims reproduce, none
mismatch, and 2 cannot be addressed from committed data.** This document records those
two, and the one reconciliation the text has already absorbed.

## 1. Resolved: the published non-failing comparison

An earlier draft of the paragraph on published non-failing comparisons reported 29 and 28 of
the 139 CS-lower genes as testable in the Chaffin DCM-vs-NF and HCM-vs-NF cardiomyocyte
tables, with 19/29 and 17/28 lower in disease against reference rates of 47% and 56%
(p=0.038, p=0.061). `scripts/spatial_first/06_published_reference.py` computes 36 and 35
testable, 25/36 and 22/35 lower in disease, reference rates 48.3% and 44.7% (p=0.009,
p=0.023). The direction of the finding is unchanged and the recomputed version is stronger:
the HCM comparison moves from marginal to conventionally significant.

The provenance of the earlier counts was never established. They were checked against every
plausible alternative definition — gene universe (all 139 lower genes; the 115 section-stable
lower core; the QC-passing subsets; the 108 lower genes testable in both platforms;
single-flag exclusions), cell type (cardiomyocyte versus pseudo-bulk rows), effect-size column
family (CellBender-corrected versus CellRanger), and the background-contamination filter on or
off. No combination yields 29 and 28 testable, and none yields 19 and 17 lower, so the
difference is not a filter choice. One value does have a likely reading: the earlier HCM
reference rate of 56% is the complement of the recomputed 44.7% lower in disease, i.e. the
fraction *not* lower, whereas the DCM rate was taken as the fraction lower.

**The current text states the recomputed values**, together with the testability caveat the
script's docstring sets out: the published tables list only genes called differential, so a
CS-lower gene absent from one of them may be unchanged in that disease or may simply not have
been reported, and the binomial reference rate is the base rate of the published differential
set rather than 0.5. These eight numbers now agree, and the harness holds the text to them.

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

## 3. Superseded material

The earlier single-nucleus arm's tables, scripts and figures are no longer in the working
tree; they are in git history at tag `pre-slim-2026-09-14`. In particular the old
`results/DE_CS_vs_*.tsv` predate the requantification the manuscript's cross-platform arm
uses and disagree with it substantially across shared genes. The tables the manuscript
relies on are `data/snrna/snrna_CS_vs_*.tsv`.

Two wordings that earlier drafts left ambiguous are now explicit in the text and checked by
the harness: the leave-one-section-out and leave-one-disease-out intersections are reported
separately (141 and 117 of the 192 hits), and the testability limit of the published tables
is stated where the 36 and 35 counts appear.
