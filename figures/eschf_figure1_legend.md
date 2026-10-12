# Figure 1. Cardiomyocyte expression in cardiac sarcoidosis by distance from immune-dominant spots and across platforms

All panels use the primary contrast: cardiac sarcoidosis (CS) against comparator cardiomyopathies over **all**
cardiomyocyte (CM)-dominant Visium spots, with no exclusion by distance to immune-dominant spots (11 CS
sections from 8 patients; 6 comparator sections from 6 patients; section-level
limma-voom; 1 CS section with no CM-dominant spot passing the UMI floor is not included). At patient level
(sections summed per patient; 8 CS and 6 comparator patients), 193 of the
223 section-level hits remained at FDR <5% and 223 of 223 kept their sign. The earlier fit restricted
to CM-dominant spots ≥0.5 mm from an immune-dominant spot is reported as a sensitivity analysis.

**(A)** H&E of one CS section (CS_108-1) and one comparator section (ARVC_1, arrhythmogenic cardiomyopathy),
with Visium spots that pass the 2000-UMI floor coloured by dominant lineage. Green, CM-dominant spots (all analysed;
1218 in CS_108-1, 2349 in ARVC_1); orange, immune-dominant spots (178 in
CS_108-1; 70 in ARVC_1). Spots of other dominant lineages are unmarked. Scale bars, 1 mm.
**(B)** Score of genes lower in CS (top; n = 116) and *ACADVL* (bottom) in CM-dominant CS spots, by distance to
the nearest immune-dominant spot. The score averages the 116 QC-passing genes lower in CS in the primary fit; its
level difference from comparators reflects how the genes were selected, and the panel shows how it varies with
distance. Score, mean log1p CP10K over member genes on depth-matched counts. Thin lines, per-patient bin means; thick
line and band, mean across patients and 95% t-interval. Blue line and band, comparator mean and range of per-patient
means. Counts beneath the axis are spots/patients per bin. Only populated bins are shown. The maximum distance from
any CS CM-dominant spot to the nearest immune-dominant spot was 1.64 mm; no CS section contained CM-dominant
myocardium beyond 2 mm from an immune-dominant spot. Text above each plot, per-patient slope: for each CS patient
with CM-dominant spots in at least two populated bins (n = 8), ordinary least squares on that patient's spots,
outcome ~ distance + CM-structural score + z(log10 UMI), with a section indicator for patients with more than one
section; shown are the mean of the per-patient slopes, its 95% t-interval and the number of negative slopes (score
+0.021 per mm, 95% CI -0.021 to +0.063; *ACADVL* +0.062 per mm,
95% CI -0.212 to +0.337). Pathway scores defined independently of the fit (KEGG TCA cycle, fatty
acid degradation and oxidative phosphorylation) were not lower in CS once fit hits were excluded.
**(C)** Spatial log2 fold change (CS vs comparators) against pooled snRNA log2 fold change (CS vs DCM and CS vs ARVC,
inverse-variance pooled) for spatial hits (FDR <5%) passing locus QC and testable in snRNA (n = 170). Because the
CS and comparator single-nucleus cohorts come from different studies, pooled single-nucleus effects beyond ±3
log2 were treated as potentially affected by cross-study quantification. Genes with pooled single-nucleus effects beyond ±3 log2 (n = 5) are shown at the axis
edge and excluded from the agreement fractions: *CFD* (-5.5), *ENG* (+4.5), *HSPB2* (+7.1), *CALM2* (+7.4), *NDUFC2* (+8.1). All five screened genes were near-absent from one
source study's quantification while present in the others, as previously observed at the FPGT-TNNI3K locus; the
±3 log2 threshold is a pre-stated sensitivity screen, not a biological limit. Triangles point outward, are drawn open in grey without a sign
category and are labelled with their value only; no other gene lies beyond the ±3 y-axis. Labelled, the genes of
(D) that pass locus QC and the screen. Text, fraction of genes with the same sign in each direction after the screen, beside the background fraction
from the label-permutation null (lower in CS 67%, n = 108, vs 48%; higher in CS 42%, n = 57,
vs 51%; without the screen 65% vs 48% and 43% vs 51%); direction agreement is reported
descriptively. ρ, Spearman correlation of fold changes after the screen.
**(D)** Concordant genes of the section-stable core; genes failing locus QC are marked. Shown
across contrasts are the 19 genes of the section-stable core (spatial FDR <5% in the full fit and in each of the
6 refits leaving out one comparator section) with the same sign and snRNA FDR <5% in the pooled single-nucleus
contrast, after the screen in (C) (*ENG* removed). Colour, direction of log2 fold change; dot size, −log10 FDR (capped at 10). Open squares, gene not
tested, or absent from the published differential-gene list (Chaffin et al. 2022, cardiomyocytes; genes with the
background-contamination flag excluded). †, locus-quantification outlier flag (*SPTAN1*); ‡, gene outside the
locus-QC universe, not evaluated by locus QC (none).

## Alt text
Four-panel figure. A: two H&E sections with coloured spots marking analysed cardiomyocyte-dominant spots and
immune-dominant spots. B: two line plots, a score of genes lower in sarcoidosis and ACADVL, against distance from the
nearest immune-dominant spot over 4 distance bins up to 2 mm, with sarcoidosis patient means and the comparator
mean and range. C: scatter of spatial against single-nucleus fold changes on a limited vertical axis, with 5
screened genes marked by open grey triangles at the edge. D: dot matrix of 19 section-stable concordant genes across five
contrasts, with many open squares in the published-reference columns.

## Graphical abstract (11 × 18 cm, width × height)
Panels B (score of genes lower in CS, *ACADVL*) and C (sign agreement after the snRNA screen) of Figure 1, redrawn with
the same populated bins and off-scale rule; off-scale points carry their value only, and *ACADVL* is the only gene
labelled.
