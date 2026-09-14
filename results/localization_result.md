# §15 Localization — does the down-in-CS program vary with distance from granuloma?

Pre-specified in `res/localization_config.json` (frozen before any fit). Per-spot scores
exported on bigpurple, job `27e6f509`, from the same raw sections, inner-join gene space
(17619 symbols) and per-section thinning (random_state=0) that produced the
published pseudobulk.

## 15.1 The stop rule fires: the positive control is flat in the retained range

| estimate | β per mm | 95% CI | p |
|---|---|---|---|
| cardiomyocyte-structural score, retained spots (≥0.5 mm) | +0.019 | -0.003, +0.041 | 0.09 |
| cardiomyocyte-structural score, **all** cardiomyocyte-dominant spots (0–1.73 mm) | +0.101 | +0.090, +0.112 | <1e-6 |

Cardiomyocyte content does rise away from lesions — strongly, and reproducing the direction
of prior work — but **only across the band the analysis excludes**. Inside the retained
range the control is flat. The 0.5 mm lesion-exclusion radius removes precisely the
distances over which composition varies, so the retained range is not a design in which a
localization gradient could be demonstrated. Per the handoff's own stop rule, this is
reported before the program result.

## 15.2 The program's pooled slope is between-section, not within-section

| estimate | β per mm | 95% CI | p |
|---|---|---|---|
| pooled distance, as pre-specified | +0.046 | +0.028, +0.064 | <1e-6 |
| distance centred on its section mean | -0.0035 | -0.022, +0.015 | 0.71 |

The pre-specified model returns a positive, highly significant slope. It does not survive
any of the checks that separate a within-tissue gradient from between-section differences:

- **Within-section distance is flat** (-0.0035 per mm, p=0.71). The pooled slope reflects
  sections whose retained spots sit farther from immune-dominant spots also having higher
  program scores overall — a between-section contrast, with 10 sections and 7 patients.
- **No section carries it.** 9 sections have ≥15 retained spots; 0 of 9 per-section slopes
  reach p<0.05 and 6 of 9 point estimates are negative (range -0.35 to +0.01 per mm).
- **The sign flips with sequencing depth.** T=1500: +0.016 (p=0.035);
  T=2000: +0.046 (p<1e-6); T=3000: -0.006 (p=0.24).
- **Dropping one section removes significance.** Leave-one-section-out: CS_102-1 p=0.12,
  CS_102-2 p=0.064; the remaining eight stay p<0.001.
- **The secondary definition reverses it.** Over all cardiomyocyte-dominant spots the program
  slope is -0.033 per mm (p<1e-6) — opposite in sign to the retained-range estimate.

The full stable-down set of 115 genes behaves identically to the 70-gene QC-passing primary
(+0.047 per mm), so the two program definitions do not differ.

## 15.3 The bound

Distance SD among retained spots = 0.199 mm; program-score SD = 0.195.

On the within-section estimate — the one not contaminated by between-section composition —
the 95% CI excludes any effect larger than **0.022 program-score units per mm**, i.e.
**0.0043 units per SD of distance = 0.022 of a program-score SD**. The effect detectable
at 80% power with this design and realised residual variance (0.116 residual SD) is
**0.026 per mm = 0.026 program-score SD per SD of distance**. A gradient of the magnitude prior
work reported for cardiomyocyte structural genes would have been detected; none is present.

## 15.4 Distance distribution in the retained set

4337 retained spots, 10 of 12 CS sections, 7 patients (7 sections with ≥50 spots).
Range 0.50–1.73 mm, median 0.69, IQR 0.27 mm. The
largest distance attained by **any** cardiomyocyte-dominant spot in the cohort is 1.73 mm.
Two sections retain nothing: CS_103 (every cardiomyocyte-dominant spot within 0.46 mm of an
immune-dominant spot) and CS_109 (no cardiomyocyte-dominant spot at all). Per-section
detail in `res/localization_persection.tsv` and `res/localization_distance_distribution.tsv`.

## 15.5 Departures from the handoff specification

1. **Depth covariate is native UMI, not thinned UMI.** The frozen spec named
   z(log10 UMI_thinned); after thinning every retained spot has exactly T counts, so that
   term is constant by construction. Pre-thinning depth is used instead.
2. **Distance is to the nearest immune-dominant spot.** No histologic granuloma annotation
   exists anywhere in this project — the region definition is compositional. The handoff's
   requested annotation-versus-proxy comparison therefore cannot be made, and the two
   distance metrics of its §5 collapse to one.
3. **The 3.97 mm tissue-scale maximum could not be reproduced or sourced.** That value does
   not appear in either repository; the observed maximum is 1.73 mm.
4. **Prior per-SD benchmarks (+0.138 structural, +0.034 program) are not unit-comparable.**
   `distance_umi_covariate.py` z-scored *pixel* distance within section, scored gene sets
   with `score_genes`, restricted to "preserved"-zone spots and applied no exclusion radius
   (33,628 spots, 8 patients). Direction is comparable; magnitude is not.
5. **Prior per-spot files named in §6 are absent** (`granulomafree_perspot.tsv.gz`,
   `distance_umi_perspot.tsv.gz`, `section_zone_summary.tsv`); only the committed
   `distance_umi_covariate.tsv` summary exists, so no prior per-spot values were re-used.
6. **The control uses 9 of 10 structural genes** — MYH6 is not in the inner-join gene space.
7. **Primary set is all retained spots (10 sections).** The project's ≥50-spot evaluable
   rule is a section-level rule; applied here as a sensitivity it gives +0.035 per mm.

## 15.6 Recommendation for paragraph 7a

The paragraph cannot be written as a clean null, because the pre-specified model returns a
significant positive slope that is an artefact of between-section composition, and the
positive control is flat in the range where the test is run. Two defensible options:

**(a) Cut it.** Preferred. The retained range (IQR 0.27 mm) does not support a
localization claim in either direction, and the flat control says so directly.

**(b) Keep it as an explicitly range-limited bound** (~90 words):

> Within cardiomyocyte-dominant spots at least 0.5 mm from immune-dominant tissue
> (4337 spots, 10 sections, 7 patients), the reduced cardiomyocyte program showed no
> gradient with distance once section identity was accounted for (β = -0.003 score units
> per mm, 95% CI -0.022 to +0.015), bounding any within-section effect below
> 0.02 SD per SD of distance. The cardiomyocyte-structural positive control was
> likewise flat over this range (p=0.09) while rising steeply across the excluded
> lesion-adjacent band (+0.10 per mm), so this range bounds — rather than tests —
> spatial localization of the program.

## 15.7 What this does not bound

- Distances below 0.5 mm; the exclusion radius was not relaxed (pre-specified).
- The up-in-CS direction, and any comparator disease — neither was tested.
- Within-section spatial correlation is not modelled. Residual nearest-neighbour lag
  correlation is median 0.07 (max 0.38) across sections, so the reported CIs are
  optimistic; the within-section null would only widen.
- Patient-level variance is estimated from 7 patients; some refits hit the variance
  boundary (statsmodels ConvergenceWarning) and the random intercept is small relative to
  residual (0.0265 vs 0.0134).
- The score is detection-dominated: median 15 of 70 program genes are detected per
  retained spot at T=2000.
