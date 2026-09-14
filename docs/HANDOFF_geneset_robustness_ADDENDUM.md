# Addendum to HANDOFF_geneset_robustness — sarcomere-content control + interpretation calibration

**Prepared for:** Claude Science (HPC). Small, cheap follow-up to the completed robustness run
(`docs/geneset_robustness_result.md`). Same pipeline, same spots, one added gene set.

---

## 0. Framing (read first — this changes how verdicts are written)

Cardiac sarcoidosis is understudied. The evidentiary bar for this paper is the **same as the field's
own foundational CS work** — Neyazi (22 patients, single cohort, descriptive) and Liu (4 CS vs 3 ICM)
— both of which report **hypothesis-generating** observations, not batch-clean multi-cohort proof.
Findings here are to be **reported as hypothesis-generating candidates with honest caveats, not
dismissed** for failing a bar that Neyazi/Liu themselves would not clear.

Concretely: a result that reproduces the *direction* of an effect in the primary CS data, is not a
scoring artifact, and is biologically coherent **meets the bar**, even if it is not disease-specific
or not significant in every reduction. The purpose of the control below is to **calibrate the wording**
of the spatial claim, not to decide whether it is publishable.

---

## 1. Why this control

The robustness run showed the anti-paracrine / preserved-enriched spatial result (analyses A + B)
reproduces for `original4`, `tnni3k`, and `mlip_pank1` but not for `objective53` or the random control.
The genes that carry it (TNNI3K, MLIP, PANK1) are **abundant cardiomyocyte structural genes**, so the
gradient could reflect either (i) a cardiomyocyte-localized CS program or (ii) generic cardiomyocyte
**content** declining where cardiomyocytes are replaced by granuloma/scar. The completed run cannot
separate these because it never scored the canonical sarcomere genes — the definitive CM-content
markers, which the manuscript explicitly *excluded* from the signature.

This addendum adds exactly that control.

## 2. New gene set (provided)

`metadata/geneset_cmstructural.txt` — canonical cardiomyocyte sarcomere / identity genes
(TNNT2, MYH7, MYH6, TTN, ACTC1, MYL2, TNNI3, TPM1, ACTN2, MYBPC3). These are the genes the paper
removed as "cardiomyocyte lineage markers"; here they are a **positive composition control**.

## 3. What to run

Score `geneset_cmstructural.txt` through the **identical** A and B logic already in
`scripts/geneset_robustness.py` (distance-to-lesion mixed model; histologic-zone gradient,
patient-level), on the same Foong CS Visium spots. Append its row to
`results/figures/geneset_robustness/summary.tsv` and the `*_ALL7.tsv` (→ ALL8) tables. No other
change.

## 4. How to read the result — and how to word the claim either way

Both outcomes are reportable; they change only the wording.

- **If canonical sarcomere genes ALSO show the preserved-enriched, away-from-lesion gradient**
  (β > 0, preserved highest): then A/B is substantially a **cardiomyocyte-content / localization**
  effect. **Word the claim as:** *"the CM-intrinsic program is expressed in, and localizes to,
  cardiomyocyte-rich myocardium and is not induced toward granulomas"* — a cardiomyocyte-localization
  statement, framed as hypothesis-generating. **Do not** upgrade it to a CS-specific cell-autonomous
  program on the strength of A/B alone, and **do not** discard it.
- **If canonical sarcomere genes do NOT show the gradient** (flat, like the random control) while the
  signature genes do: then the signature's spatial gradient is **not** explained by CM content, and
  the stronger wording — *"a cardiomyocyte-localized program whose spatial distribution is not a
  simple consequence of cardiomyocyte abundance"* — is supported, still as a hypothesis-generating
  observation.

Either way the spatial result stays in the paper; the control decides only whether it is described as
"cardiomyocyte-localized" (weaker) or "beyond cardiomyocyte content" (stronger).

## 5. Interpretation calibration for the whole robustness table (apply to the manuscript)

Re-state the completed-run verdicts at the Neyazi/Liu bar — as candidate/hypothesis-generating, not
pass/fail dismissals:

- **A/B (spatial localization):** a reproducible, method-validated cardiomyocyte-localized,
  non-paracrine pattern. Report it; word per §4.
- **C (stage-independence):** consistent with a program present across remodeling stages, but the
  snRNA test is underpowered (n=22, control also flat) — present qualitatively, do not overstate.
- **D (cross-disease):** **hypothesis-generating.** The four-gene signature ranks CS highest across
  the available comparators; this is a candidate diagnostic direction, explicitly limited by mimic N
  (2/disease), no HCM arm, and spot-level pseudoreplication. Frame as a lead to test, not a
  validated specificity claim — the same status Neyazi/Liu give their own single-cohort findings.
- **Individual genes:** TNNI3K and GJB7 are **candidate CS-associated genes** (GJB7 was independently
  noted by Neyazi); MLIP/PANK1 are a shared arrhythmogenic-remodeling component (batch-clean vs
  ischemic in Liu). All hypothesis-generating.

## 6. Deliverable

Append the sarcomere-control row to `summary.tsv` / `*_ALL8.tsv`, add one paragraph to
`docs/geneset_robustness_result.md` under a "Sarcomere-content control" heading stating which §4
outcome held and therefore which wording the manuscript should use. Keep everything
hypothesis-generating in tone.
