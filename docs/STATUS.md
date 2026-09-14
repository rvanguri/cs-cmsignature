# Project status — CS cardiomyocyte signature

> **Superseded as the project's headline.** Everything below describes the earlier
> single-nucleus arm and its four-gene candidate panel (GJB7, TNNI3K, MLIP, PANK1). The
> current manuscript is the spatial-first brief report *"Metabolic gene expression differs
> in non-lesional myocardium between cardiac sarcoidosis and other cardiomyopathies"* —
> see `README.md` for its result and `scripts/spatial_first/run_all.sh` for the pipeline
> that regenerates it. This file is kept because its limitations section still governs how
> the single-nucleus data may be used, and because the cross-platform arm of the current
> manuscript draws on the same cohort.

**Updated:** 2026-09-04. One-page state for whoever picks it up next.

**Evidentiary bar.** Cardiac sarcoidosis is understudied; the foundational CS transcriptomic papers
(Neyazi, 22 patients, single cohort; Liu, 4 CS vs 3 ICM) report **hypothesis-generating**, honestly
caveated observations, not batch-clean multi-cohort proof. This project is held to the **same bar**:
findings are reported as candidates with clear limitations, not dismissed for failing a standard the
field itself does not meet.

**TL;DR.** The paper's core is a **candidate four-gene cardiomyocyte-associated signature of CS**
(GJB7, TNNI3K, MLIP, PANK1), hypothesis-generating. Two refinements the analysis established:
MLIP/PANK1 behave as a *shared* arrhythmogenic-remodeling component (not CS-exclusive), and the
spatial result is best worded as **cardiomyocyte-localized / not paracrine-induced** rather than
"cell-autonomous." TNNI3K and GJB7 remain CS-leaning candidates (GJB7 independently noted by Neyazi).

---

## Established (facts, model-free or batch-clean)

- **MLIP and PANK1 are shared with arrhythmogenic cardiomyopathy.** In the raw pairwise snRNA
  comparison they are equal to (MLIP, p=0.83) or below (PANK1) DCM/ARVC; batch-clean support is vs
  *ischemic* cardiomyopathy only (Liu: MLIP padj 0.045, PANK1 padj 0.019). Best described as a shared
  arrhythmogenic-remodeling component, not a CS-exclusive marker.
- **The four-gene panel was curated, not top-ranked.** By the pipeline's `evidence_score` the four
  rank 336 / 420 / 1,978 / 3,781 of ~18,700 (CTNNA3, a canonical ARVC gene, is #1). Selection used a
  biological "relevance to arrhythmia" filter — legitimate, but the panel should be presented as an
  illustrative candidate set, not an optimized or unique one.
- **The cross-cohort discovery contrast is confounded** (rank-deficient; housekeeping genes score
  "CS-specific" even after quantile normalization). Cross-cohort DE is therefore hypothesis-generating,
  as the manuscript already states — not a specificity proof.
- **A batch-clean-concordant signature exists and is the granuloma/fibrosis program.** The 111-gene
  `high_confidence_cs_set` (significant + concordant across every primary contrast **and** Liu) is
  dominated by immune/myeloid (CD163, CD74, CSF2RA…) and collagen/ECM (COL1A1/3A1, FN1) genes — the
  known CS inflammatory/fibrotic biology. The cardiomyocyte panel is the residual after removing those
  lineage genes. Disease-distance clusters CS with DCM/ARVC (~12–13) and isolates HCM (74).
- **Provenance.** Neyazi CS is stage-heterogeneous (VAD n=7, transplant n=13, autopsy/myectomy n=2;
  EF≥35 in 4/22); all Visium (CS + mimics) is end-stage transplant.

## Candidate / hypothesis-generating (report with caveats)

- **Candidate four-gene CS cardiomyocyte signature** (GJB7, TNNI3K, MLIP, PANK1). Nominated from the
  (confounded) discovery + Foong replication + curation. Illustrative, not unique.
- **Spatial: cardiomyocyte-localized, not paracrine-induced.** The program is expressed in
  preserved, cardiomyocyte-rich myocardium and is not enriched toward granulomas (distance β=+0.044,
  8/8 patients). **Gene-set robustness (ran) shows this is a cardiomyocyte-localization/content
  effect:** canonical sarcomere genes (TNNT2, MYH7, TTN…) show the *same, stronger* gradient (β=+0.162,
  8/8; preserved≫fibrotic≫granuloma, p=0.004). Method validated — inflammation/fibrosis controls score
  the opposite polarity (`highconf111` β=−0.044, a near-mirror), negative control fails. **Word it as
  "cardiomyocyte-localized / not paracrine," not "cell-autonomous CS-specific."**
- **Stage-independence within CS.** Consistent (flat across zone/fibrotic-burden/procedure/EF), but
  the snRNA test is underpowered (n=22; random control also flat) — present qualitatively.
- **Cross-disease direction (candidate diagnostic lead).** The four-gene signature ranks CS highest
  across available comparators (CS > control > ARVC > LMNA), but this depends on the exact panel
  (carried by GJB7; components rank CS lower) and is limited by mimic N (2/disease), no HCM arm, and
  spot-level pseudoreplication. A lead to test, not a validated specificity claim — the same status
  Neyazi/Liu give their single-cohort findings.
- **TNNI3K, GJB7 — CS-leaning candidate genes.** GJB7 was independently reported upregulated in CS
  cardiomyocytes by Neyazi. TNNI3K has the largest CS-vs-mimic effect but no batch-clean confirmation
  (see Open).

## Open (genuinely undecided)

- **TNNI3K: CS-associated or length/batch artifact?** 64× higher in CS than ARVC in snRNA (does not
  fit "shared"), but unconfirmed batch-clean and a 309 kb intron-rich gene prone to intronic-counting
  inflation. Pivotal; see `REFRAME_honest_fallback.md` §5.1. Resolved only by matched-pipeline
  re-quantification (blocked) or an independent CS cohort.
- **Does the CM signature do anything beyond cardiomyocyte content spatially?** The sarcomere control
  says the localization gradient is largely CM content; whether the signature adds CS-specific spatial
  information beyond that is not established.

## Blocked

- **Unified re-quantification** (would adjudicate TNNI3K and cross-cohort specificity): needs raw
  Reichart FASTQ from EGA `EGAS00001006374`, access-controlled and not obtainable.
  `HANDOFF_unified_requantification.md` is parked.

## Done this cycle

- **Gene-set robustness check — RAN** (`docs/geneset_robustness_result.md`, `results/figures/
  geneset_robustness/`), incl. the sarcomere-content control addendum. Conclusion: spatial A/B is a
  validated cardiomyocyte-localization effect; C underpowered; D panel-dependent/hypothesis-generating.
- `CS_cardiomyocyte_JACC_BriefReport_v2.docx` — reframe (arrhythmogenic-class). *Being superseded by v3
  with the cardiomyocyte-localized / hypothesis-generating wording.*
- `docs/REFRAME_honest_fallback.md` (incl. §5.1 TNNI3K, §5.2), `docs/HANDOFF_geneset_robustness.md`
  (+ addendum, + gene-set files), `docs/HANDOFF_unified_requantification.md` (blocked, for the record).

## Recommended framing for the manuscript

A hypothesis-generating Brief Report: a candidate cardiomyocyte-associated CS signature that is
**cardiomyocyte-localized and not paracrine-induced**, sits within a shared arrhythmogenic-remodeling
program (MLIP/PANK1 shared; TNNI3K/GJB7 CS-leaning candidates), and ranks CS highest across a small
disease panel as a candidate diagnostic lead — all limited by single-cohort discovery and small mimic
N, and offered in the same hypothesis-generating spirit as the field's foundational CS studies.
