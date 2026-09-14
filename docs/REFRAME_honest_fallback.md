# Reframing plan — the honest fallback (stands on batch-clean evidence, no re-quantification)

**Companion to** `docs/HANDOFF_unified_requantification.md` (the blocked *upgrade* path) and
`docs/geneset_robustness_result.md` (the robustness check, now run). **This** is what the paper can
defensibly claim right now.

**Evidentiary bar.** Cardiac sarcoidosis is understudied; Neyazi (22 patients, single cohort) and Liu
(4 vs 3) report **hypothesis-generating** findings, not batch-clean multi-cohort proof. This reframe
holds the paper to the **same bar** — it narrows over-specific claims and words them precisely, but it
**reports candidates rather than dismissing them.** Where a claim below is called "not supported,"
that means *not established as specific/mechanistic*, not *discard the observation*.

---

## 1. The core problem this reframing solves

The current manuscript claims a **cell-autonomous four-gene signature (GJB7, TNNI3K, MLIP, PANK1)
that is CS-specific and points toward molecular diagnosis.** Re-analysis shows the specificity claim
does not hold:

- The cross-cohort discovery contrast is confounded at the gene level (CS is a single cohort;
  ribosomal/OXPHOS housekeeping genes score as "CS-specific" even after quantile normalization).
- **MLIP and PANK1 are not CS-specific** vs the arrhythmogenic mimics (raw CS≈DCM, p=0.83; CS lower
  than ARVC); they separate CS only from HCM/NF.
- **TNNI3K and GJB7** survive as candidates but have no batch-clean confirmation.

So the honest paper stops claiming specificity/diagnosis and instead reports what the **batch-clean,
within-cohort** evidence actually supports.

---

## 2. New thesis

> In cardiac sarcoidosis we nominate a candidate four-gene cardiomyocyte-associated program
> (GJB7, TNNI3K, MLIP, PANK1) that is **cardiomyocyte-localized and not paracrine-induced** and is
> present across remodeling stages. Across diseases it behaves as part of an
> **arrhythmogenic-cardiomyopathy remodeling program** (shared with ARVC and LMNA cardiomyopathy,
> enriched over hypertrophic cardiomyopathy and normal myocardium); MLIP and PANK1 are its shared
> component, while TNNI3K and GJB7 are CS-leaning candidate genes (GJB7 independently noted by Neyazi)
> that warrant further study.

It is a **hypothesis-generating characterization-and-nomination** paper — pitched at the same bar as
the field's foundational single-cohort CS studies — not a validated specificity/diagnostic paper.

---

## 3. What the spatial (Foong) arm actually shows — the reframed reading

The Foong series (GSE314910) is **CS-only, 12 samples / 9 patients**, end-stage. Two spatial results,
correctly interpreted:

**(a) Cardiomyocyte-localized, not paracrine-induced.** Within-patient, the CM program rises with
distance *away* from granulomas (mixed model β=+0.044/SD, CI +0.041 to +0.046; positive slope in 8/8
patients) — it peaks in preserved cardiomyocyte tissue, not at lesions, so it is not a paracrine
granuloma response. The gene-set robustness run (now complete) validated the method (inflammation/
fibrosis controls score the opposite polarity; random control fails) **and** showed this gradient is
largely a **cardiomyocyte-localization/content** effect: canonical sarcomere genes (TNNT2, MYH7, TTN…)
show the same, stronger gradient (β=+0.162). So word this as *cardiomyocyte-localized and
non-paracrine* — a genuine, hypothesis-generating observation — **not** as a CS-specific
cell-autonomous program. It does not depend on the cross-cohort contrast.

**(b) Genuine in-situ expression, not a dissociation artifact.** TNNI3K and GJB7 are truly expressed
in CS cardiomyocyte-dominant spots on an independent platform, rebutting the concern that a
nuclear-retained long gene (TNNI3K, 309 kb) or GJB7 are single-nucleus capture artifacts.

**(c) The cross-disease pattern is the informative part.** Per-gene, CM-dominant spots:

```
CS  ≈  ARVC  ≈  LMNA   >   HCM  ≈  normal
```

- Enriched over **Kuppe normal** → disease-associated remodeling genes (TNNI3K 0.52 vs 0.19;
  GJB7 0.009 vs 0.000; MLIP 1.16 vs 0.38; PANK1 0.135 vs 0.025).
- Enriched over **HCM** (Foong-published CS-vs-HCM Visium DE: GJB7 log2FC 5.96, padj 3.6e-9; TNNI3K
  0.75, padj 7.5e-4).
- **Not** enriched over **ARVC / LMNA** (TNNI3K: CS 0.52 vs DC-LMNA 0.81, ARVC_2 0.50; GJB7: CS 0.009
  ≈ ARVC_2 0.009).

**Interpretation.** This shape rules out "disease vs health" (HCM would be up) and "cardiomyopathy in
general" (HCM stays at baseline). What CS, ARVC, and LMNA share and HCM/normal lack is
**arrhythmogenic / dilated-type cardiomyocyte remodeling** — the structural–electrical remodeling axis
of the diseases that cause ventricular arrhythmia, conduction disease, and sudden death. HCM is a
hypertrophic/diastolic disease with different remodeling biology, so it sits with normal. The gene
identities fit exactly: **GJB7** (gap junction → electrical/conduction remodeling), **MLIP** (LMNA
pathway → shared with LMNA cardiomyopathy by mechanism), **TNNI3K** (sarcomeric/structural
remodeling), **PANK1** (metabolic remodeling).

**So the signature marks the arrhythmogenic-cardiomyopathy class; CS cardiomyocytes molecularly
resemble those of ARVC and LMNA and are distinct from HCM and normal — it is not CS-specific.**

> Do **not** call the Kuppe enrichment "late-stage remodeling." It conflicts with our own batch-clean
> results: within CS the program is stage-independent (flat across fibrotic burden, EF, procedure)
> and, spatially, is **highest in preserved zones and lowest in fibrotic/granulomatous** — the
> opposite of a scar/severity marker. Describe it as a disease-associated program present across CS
> stages, including preserved tissue.

---

## 4. Evidence tiers (how to present)

**Tier 1 — batch-clean / within-cohort (the spine):**
- Spatial cell-autonomy (Foong distance model) — cardiomyocyte-intrinsic, anti-paracrine.
- Stage-independence within CS (Neyazi snRNA zones + patient stage; Foong preserved-vs-scar).
- Liu CS-vs-ICM within-study: MLIP padj 0.045, PANK1 padj 0.019 (CS > ischemic).

**Tier 2 — hypothesis-generating / confounded (labeled as such):**
- Cross-cohort CS-vs-DCM/ARVC/HCM/NF → nominates TNNI3K/GJB7 candidates; confound disclosed
  (single CS cohort; quantification batch; housekeeping negative-control failure).

**Descriptive / small-series:**
- The n=1–2 arrhythmogenic Visium mimics (ARVC, LMNA) → establish the arrhythmogenic-class pattern
  in §3(c) as **directional**, and motivate the unified re-quantification; not a powered test.

---

## 5. Gene-by-gene honest status

| Gene | Batch-clean support | Cross-cohort | Honest label |
|---|---|---|---|
| TNNI3K | none testable (Liu 5′ under-captures) | large CS>all mimics (confounded) | **unresolved — see §5.1**; either the one real CS-specific gene or a length/batch artifact |
| GJB7 | none testable (floor) | modest CS>all mimics; strong Foong CS-vs-HCM | CS-associated **candidate**; arrhythmogenic-class in spatial |
| MLIP | **CS > ICM (Liu, p=0.045)** | equal to DCM/ARVC | **shared** arrhythmogenic-remodeling (LMNA pathway); not CS-specific |
| PANK1 | **CS > ICM (Liu, p=0.019)** | lower than DCM/ARVC | **shared** arrhythmogenic-remodeling; not CS-specific |

Tension to state plainly: the gene with the strongest CS signal (TNNI3K) has the weakest
confirmation; the genes with batch-clean confirmation (MLIP/PANK1) are not CS-specific. Hence the
paper leads on **compartment-level** biology and treats genes as candidates.

---

## 5.1 Pivotal gene: TNNI3K (unresolved — do not resolve in either direction)

TNNI3K breaks the shared-arrhythmogenic pattern that MLIP/PANK1 fit. A gene shared across
arrhythmogenic cardiomyopathies should be roughly *equal* in CS versus DCM/ARVC in the snRNA; TNNI3K
is instead ~64× higher in CS than in ARVC (raw MWU p<1e-8). So it does **not** behave as shared in
the snRNA — which forces a fork with two mutually exclusive readings, and current data cannot
distinguish them:

- **A — genuinely CS-specific.** The snRNA elevation is real; TNNI3K is more up in CS than in any
  other cardiomyopathy, including arrhythmogenic ones. Then it does **not** belong in the "shared"
  bucket, and its "shared" appearance in the spatial data is because a 309 kb gene is barely detected
  by poly-A Visium (it floors out and cannot show the difference). Under A, TNNI3K is the project's
  one real CS-specific finding.
- **B — a length/batch artifact.** TNNI3K is long and intron-rich, and Neyazi quantified with
  intronic reads included (per its GEO record); long intron-rich genes inflate under
  intronic-inclusive quantification, so a Neyazi-vs-comparator pipeline difference could manufacture
  the 64×. The same 309 kb length explains why poly-A Visium cannot see it. Under B, TNNI3K is
  arrhythmogenic-shared like MLIP/PANK1 and nothing special.

Evidence splits: toward **B** — it does not replicate in the batch-clean Liu contrast (capture-limited
but null), its length makes it the prime intron-count suspect, and a 64× snRNA gap collapsing to
*equal* in spatial is larger than detection-compression alone comfortably explains; toward **A** — the
effect is enormous with non-overlapping distributions in a cardiac-specific kinase, and MLIP is
*equally long* (338 kb) yet is **not** CS-inflated, arguing the batch is not purely length-driven and
leaving room for real signal.

**Instruction for the manuscript:** do neither tempting thing — do not fold TNNI3K into
"shared arrhythmogenic" (the snRNA contradicts that) and do not claim it as CS-specific (unconfirmable;
batch plausible). Flag it explicitly as the pivotal unresolved gene. **What resolves it:** unified
re-quantification of CS + DCM/ARVC through one pipeline with matched intronic handling (see the
re-quantification handoff — now blocked by EGA access), an independent CS cohort, or, partially, the
TNNI3K-only arm of the gene-set robustness check (interpreting a null there against the poly-A
detection caveat).

## 5.2 Retired defense: "the snRNA DE supports CM-intrinsic specificity"

`HANDOFF_pergene_crossdisease-2.md` §9 defends the signature by conceding the spatial CM-intrinsic
score is at baseline (§9.5) and falling back on the snRNA arm: specificity "rests on the snRNA
pseudobulk arm — CS vs DCM/ARVC/HCM/NF, 142 units / 122 patients … where the CM-intrinsic genes are
strongly DE" (§9.3). **This fallback is retired; do not use it anywhere in the manuscript.**

- The snRNA contrast is rank-deficient/confounded (the DE table flags `contrast_confounded=True`), and
  ribosomal/OXPHOS housekeeping genes score as "CS-specific" **even after quantile normalization** —
  so "strongly DE in the snRNA arm" is not evidence of CS specificity for anything.
- In the model-free raw pairwise comparison within that same snRNA, MLIP and PANK1 are equal to (or
  lower than) DCM/ARVC — not CS-specific there either.
- The objective CS-specific set (led by TNNI3K) is **disjoint** from the batch-clean set: 0 of 53
  genes replicate in Liu.

The handoff is also internally inconsistent: §9.2 calls the CM-intrinsic axis "*by construction the
shared structural-remodeling axis*," which cannot coexist with §9.3's claim that its between-disease
specificity rests on the snRNA DE. We resolve the contradiction toward **shared**.

Where we *agree* with that handoff, and carry forward: the CM-intrinsic axis is shared with
arrhythmogenic cardiomyopathies (§9.2); the spatial CM-intrinsic numbers are baseline and must not be
presented as positive evidence (§9.5); spatial CS-specificity lives in the inflammation genes
(NLRC4/IL7), which are a curated positive control. The net effect is only to remove the one escape
hatch — the snRNA DE cannot be used to keep a CS-specific CM claim alive.

---

## 6. Edits to make (drop / keep / add)

**Drop:** "molecular diagnosis," "diagnostic assay," "biopsy yield," and any claim of established
**CS-specificity**.

**Reclassify:** the four genes into "shared arrhythmogenic-remodeling (MLIP/PANK1)" vs
"CS-associated candidate (TNNI3K/GJB7)"; recast the whole panel as an **arrhythmogenic-cardiomyopathy
cardiomyocyte program**, not a CS signature.

**Demote:** the cross-cohort log2FC figure/claim to hypothesis-generating, with the confound in the
caption.

**Add:** the stage-independence result; the arrhythmogenic-class spatial interpretation (§3c); the
explicit confound disclosure (single CS cohort + housekeeping negative controls currently DE); and
the cross-disease positioning (CS resembles ARVC/LMNA, distinct from HCM/normal).

**Keep unchanged:** cell-autonomy (Foong distance model) and the within-CS spatial/stage
characterization — never depended on the cross-cohort contrast.

---

## 7. Figure implications

Lead with the batch-clean within-CS spatial story (preserved-vs-scar + distance/cell-autonomy) and
the cross-disease pattern (CS ≈ ARVC/LMNA > HCM ≈ normal) with n annotated. Move the cross-cohort
heatmap to a supporting/hypothesis-generating role with a confound caption. No claim panel a reviewer
can overturn by running the pairwise test.

---

## 8. Residual caveats (carry in Limitations)

- TNNI3K/GJB7 have **no batch-clean confirmation** (candidates only); TNNI3K is poly-A
  detection-compressed and GJB7 sits at the spatial floor.
- The arrhythmogenic-class comparison rests on **n=1–2** ARVC/LMNA Visium sections → directional.
- HCM being low may reflect less-remodeled tissue (myectomy/preserved EF) as well as disease type;
  the ARVC/LMNA equivalence is what tips interpretation toward disease-class over severity.
- Stage-independence shown via within-heart histologic zones from already-advanced hearts, not
  true early-timepoint biopsies.
- The inflammation axis (NLRC4/IL1RAP/BACH2/IL7) is a **curated positive control**, not a discovery.

---

## 9. Venue consequence

Less "diagnostic punch," more mechanistic/taxonomic: CS cardiomyocyte state in situ + its position
within the arrhythmogenic-cardiomyopathy spectrum + candidate nomination. Fits a mechanistic or
spatial-cardiac-biology venue better than a diagnosis-framed one. Pairs with the re-quantification
handoff: this is what stands now; the re-quant is what could later upgrade TNNI3K/GJB7 from candidate
to established.
