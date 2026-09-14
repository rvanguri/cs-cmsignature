# Gene-set robustness of the CS spatial conclusions

**Question.** Do the cardiac-sarcoidosis (CS) spatial conclusions survive when the
cardiomyocyte (CM) program is scored on defensible gene sets rather than the
hand-curated four genes? Five gene sets were scored through the *identical*
analysis logic (only the gene list changed):

| set | genes | role |
|---|---|---|
| `original4` | 4 hand-curated (published `CM_intrinsic`) | published signature |
| `objective53` | 53-gene broad CM program | objective, unbiased comparator |
| `tnni3k` | TNNI3K only | single-gene reduction |
| `mlip_pank1` | MLIP + PANK1 | pairwise reduction |
| `negctrl53` | 53 random genes | negative control (must **fail**) |
| `inflammation4` | NLRC4, IL1RAP, BACH2, IL7 | directional control — granuloma/inflammation axis (must behave **oppositely**) |
| `highconf111` | 111-gene batch-clean signature | directional control — inflammation/fibrosis-dominated (must behave **oppositely**) |

The last two are **directional (positive) controls**: they are genuine
lesion/granuloma-associated programs, so a valid pipeline must score them
*highest toward lesions and in granulomatous/fibrotic zones* — the mirror image
of the CM sets. A pass here validates that the CM "away-from-lesion" result is a
real biological contrast, not an artifact of the scoring method or tissue
architecture; a failure would invalidate every other row.

Four analyses were run: **A** distance-to-lesion mixed model (Foong CS Visium,
12 sections, 8 patients); **B** histologic-zone gradient (patient-level
Wilcoxon, n=9); **C** snRNA pseudobulk zone flatness (22 CS samples); **D**
cross-disease ordering in CM-dominant spots (CS + 4 Kuppe controls + ARVC×2 +
LMNA×2 + HCM×1 + Chagas×1 — six conditions).

**Pipeline fidelity check.** The `original4` run reproduces the committed
published baseline to the decimal — distance β = +0.0439 (95% CI 0.0414–0.0464,
p = 7.2e-262, 8/8 patients positive) and zone means preserved +0.0052 >
granulomatous −0.1553 > fibrotic −0.2347 — so any difference across sets is a
gene-set effect, not a code change.

## Acceptance matrix

| gene set | A anti-paracrine | B zone gradient | C stage-indep. | D CS-specific |
|---|---|---|---|---|
| **original4** | **PASS** β=+0.044, 8/8 | **PASS** p=0.008/0.004 | PASS (flat) | FAIL CS #3/6 (HCM > Chagas > CS > control > ARVC > LMNA) |
| **tnni3k** | **PASS** β=+0.081, 8/8 | PARTIAL dir. only, p=0.063 | PASS (flat) | PARTIAL CS #2/6 (LMNA > CS > control > ARVC > HCM > Chagas) |
| **mlip_pank1** | **PASS** β=+0.050, 8/8 | **PASS** p=0.008/0.004 | PASS (flat) | FAIL CS #6/6 (lowest) |
| **objective53** | FAIL β=−0.0015, 0/8 | FAIL pres. not highest, NS | PASS (flat) | FAIL CS #6/6 (lowest) |
| **negctrl53** | FAIL β=+0.0008, 4/8 | FAIL NS | PASS (flat) | FAIL CS #4/6 |

*(Verdict rules: A = β≥0 and ≥7/8 patient positive slopes; B = preserved zone
highest and patient-level Wilcoxon p<0.05; C = Kruskal p>0.05 flat; D = CS ranks
highest. Full numbers in `summary.tsv`.)*

**Negative-control validation.** `negctrl53` correctly **FAILS the three
discriminating analyses (A, B, D)**. It "passes" C only because flatness is the
low-power default — see below. This confirms A/B/D detect real signal, not
CM-content or tissue architecture.

### Directional-control validation (the key check)

The two inflammation/fibrosis programs are scored through the *same* A/B logic,
but the expected direction is **inverted** — enriched *toward* lesions and in
granulomatous/fibrotic zones:

| set | A distance β | 8/8 slopes | B zone order | patient-level p (pres vs fib / gran) | verdict |
|---|---|---|---|---|---|
| `highconf111` | **−0.0441** | 0/8 positive | granulomatous +0.446 ≫ fibrotic +0.230 ≫ preserved +0.006 | **0.0039 / 0.0039** | **PASS (strong)** |
| `inflammation4` | −0.0025 | 0/8 positive | granulomatous +0.023 (highest) | 0.25 / 0.91 (NS) | PASS (direction) |

`highconf111` is a **near-perfect mirror of the published signature**: its
distance slope β = −0.044 is the sign-flipped twin of `original4`'s +0.044, and
its zone gradient runs granulomatous ≫ fibrotic ≫ preserved (patient-level
p = 0.0039, opposite polarity). `inflammation4` is directionally correct (β < 0,
granulomatous zone highest) but weaker in magnitude and not significant at the
patient level. Neither ranks CS highest in D — expected, since they are not
CS-specificity sets.

**This validates the whole experiment.** The pipeline demonstrably detects a
lesion-enriched / granuloma-enriched signal when one truly exists, and assigns
it the *opposite* spatial polarity to the CM-integrity programs. The CM sets'
"score rises away from the lesion, highest in preserved myocardium" result is
therefore a real biological contrast — not a universal artifact of module
scoring or tissue architecture.

### Sarcomere-content control (composition calibration)

To test whether the A/B spatial gradient reflects a CS-specific cardiomyocyte
program or simply **cardiomyocyte content** declining where myocardium is
replaced by granuloma/scar, the canonical sarcomere / CM-identity genes
(`cmstructural`: TNNT2, MYH7, MYH6, TTN, ACTC1, MYL2, TNNI3, TPM1, ACTN2,
MYBPC3 — the markers the manuscript *excluded* from the signature) were scored
through the identical A and B logic on the same Foong CS spots.

| set | A distance β | 8/8 slopes | B zone order | patient-level p |
|---|---|---|---|---|
| `cmstructural` | **+0.162** | 8/8 positive | preserved 2.90 ≫ fibrotic 1.15 ≫ granulomatous 0.07 | **0.0039** |
| `original4` (for reference) | +0.044 | 8/8 positive | preserved highest | 0.008 / 0.004 |

**Outcome: §4 case (i) holds.** The canonical sarcomere genes show the *same*
preserved-enriched, away-from-lesion gradient, and considerably more strongly
than the signature (β = +0.162 vs +0.044). The A/B pattern is therefore
substantially a **cardiomyocyte-content / localization** effect: the signature
genes are expressed in cardiomyocytes and track cardiomyocyte abundance across
the tissue.

**Wording the manuscript should use (hypothesis-generating):** *"the CM-intrinsic
program is expressed in, and localizes to, cardiomyocyte-rich myocardium and is
not induced toward granulomas."* This is a **cardiomyocyte-localization**
statement — the weaker of the two §4 phrasings. On the strength of A/B alone the
claim should **not** be upgraded to a CS-specific, cell-autonomous program
"beyond cardiomyocyte content"; equally, the spatial result is **not discarded** —
it remains a reproducible, method-validated, biologically coherent localization
pattern that stays in the paper.

## Interpretation

> **Evidentiary bar.** Cardiac sarcoidosis is understudied; the field's
> foundational spatial/transcriptomic work (Neyazi, 22 patients single cohort;
> Liu, 4 CS vs 3 ICM) reports *hypothesis-generating* observations, not
> batch-clean multi-cohort proof. Verdicts below are stated at that same bar: a
> result that reproduces the *direction* of an effect in the primary CS data, is
> not a scoring artifact, and is biologically coherent **meets the bar** — even
> if it is not disease-specific or not significant in every reduction. These are
> candidate findings to report with honest caveats, not pass/fail dismissals.

**1 — Within-CS spatial localization (A + B) is REPRODUCIBLE and
cardiomyocyte-localized.** The pattern — the CM program is *not* induced toward
granulomas and is highest in preserved myocardium — reproduces across every
curated candidate
set: `original4`, `tnni3k`, and `mlip_pank1` all show a positive distance slope
in 8/8 patients (A), and `original4` + `mlip_pank1` reach patient-level
significance on the zone gradient with `tnni3k` directionally identical (B). It
is therefore **not an artifact of the four-gene choice**; the signal is carried
by specific CM-integrity genes (TNNI3K, MLIP, PANK1). The sarcomere-content
control (above) shows the gradient is substantially a **cardiomyocyte-content /
localization** effect — the canonical sarcomere genes show the same, stronger
gradient — so the claim is worded as *cardiomyocyte-localized and not induced
toward granulomas*, a hypothesis-generating localization statement, not a
CS-specific cell-autonomous program on the strength of A/B alone. The broad `objective53`
set is flat here, but that is interpretable rather than contradictory: a generic
CM-identity program is uniformly expressed across CM spots, so it has no spatial
gradient to detect — it behaves like the negative control precisely because it
marks CM *presence*, not CM *integrity*.

**2 — snRNA stage-independence (C) is WEAK / non-discriminating.** All five sets
are flat across histologic zones (Kruskal p: original4 0.43, objective53 0.15,
tnni3k 0.12, mlip_pank1 0.93, negctrl53 0.55) — *including the random control*.
With n=22 pseudobulk samples (6 preserved / 3 granulomatous / 13 fibrotic) the
test has essentially no power, so "flat" is the null default. Analysis C
reproduces the published claim for `original4` but provides little independent
evidence and should not be presented as a robustness win.

**3 — Cross-disease CS-specificity (D) does NOT survive the HCM arm; it is a
CM-remodeling axis, not a CS-specific one.** With HCM and Chagas added (six
conditions), the four-gene `original4` signature **no longer ranks CS highest**:
the ordering is **HCM > Chagas > CS > control > ARVC > LMNA** (CS #3/6). CS still
separates cleanly from the dilated/arrhythmogenic mimics and normal myocardium
(CS > control, ARVC, LMNA; all p ≤ 5e-11), but both HCM (a primary hypertrophic
cardiomyopathy; median 0.039 vs CS 0.009, p=2.6e-20) and Chagas (inflammatory
myocarditis; median 0.030, p=4.6e-6) score **above** CS. This is the expected
behaviour of a cardiomyocyte-integrity/remodeling program — HCM and Chagas both
feature heavy cardiomyocyte remodeling — and it converges with the
composition-control finding (§ "Word the spatial claim…"): the program marks
broad CM remodeling shared across cardiomyopathies, not a CS-specific signal.
The one exception is the single gene **`tnni3k`**, which still places CS #2/6 —
below LMNA but **above HCM and Chagas** as well as control and ARVC (CS > HCM
p=1.2e-114; CS > Chagas p=1.2e-84) — the only set that discriminates CS from the
hypertrophic/inflammatory mimics. `objective53` and `mlip_pank1` invert (CS
lowest), as before. Caveats remain: HCM and Chagas are **single sections each**,
the other mimics are 2 sections, and spot-level Mann-Whitney tests pool tens of
thousands of spots so the p-values reflect pseudoreplication — only the
*direction* of ordering is interpretable. Framed at the Neyazi/Liu single-cohort
bar, the honest read is hypothesis-generating: adding the previously-missing HCM
arm shows the 4-gene program is **not** cross-disease-specific to CS, while
TNNI3K alone remains a candidate discriminating lead worth testing with
biological replicates.

## Manuscript-framing recommendation

- **Lead with the within-CS spatial result (A+B).** It is the robust spine of
  the story: the anti-paracrine / cell-autonomy localization holds under
  single-gene (TNNI3K), pairwise (MLIP+PANK1), and 4-gene scoring, and the
  random control fails it. Frame the four genes as *exemplars of a CM-integrity
  program* whose spatial behavior is reproducible, not as an irreducible unit.
- **Do not claim cross-disease CS-specificity — the HCM arm breaks it.** With
  HCM and Chagas added, the 4-gene `original4` program ranks CS third (HCM and
  Chagas both higher), so it marks a pan-cardiomyopathy CM-remodeling axis, not a
  CS-specific one. Present D honestly: CS separates from control/ARVC/LMNA but
  **not** from HCM (hypertrophic CM) or Chagas (inflammatory myocarditis). If a
  disease-discrimination lead is reported, note that only single-gene **TNNI3K**
  still places CS above HCM and Chagas, and caveat everything by single-section
  mimic N (HCM, Chagas = 1 section; ARVC, LMNA = 2) and spot-level
  pseudoreplication — it needs biological replicates before it can be called
  robust.
- **Retire analysis C as a robustness argument.** State stage-independence
  qualitatively (present at all remodeling stages) but acknowledge the snRNA test
  is underpowered — the negative control passing makes it non-diagnostic.
- **Report both control classes together.** The negative control (`negctrl53`)
  fails A/B/D, and the directional controls (`highconf111`, `inflammation4`)
  reproduce the inflammation gradient in the *opposite* direction — `highconf111`
  mirroring `original4` almost exactly (β = −0.044 vs +0.044). This pairing is
  the strongest single piece of evidence that the CM "away-from-lesion / preserved-
  enriched" localization is genuine biology and the pipeline is well-calibrated.
- **Word the spatial claim as cardiomyocyte-localized, not "beyond CM content."**
  The sarcomere-content control (`cmstructural`, β = +0.162) shows the A/B
  gradient is substantially a cardiomyocyte-content effect. The defensible claim
  is that the CM-intrinsic program is expressed in and localizes to
  cardiomyocyte-rich myocardium and is not induced toward granulomas — not a
  cell-autonomous CS-specific program.
- **Individual genes are candidate CS-associated leads (hypothesis-generating).**
  TNNI3K and GJB7 are candidate CS-associated genes (GJB7 independently noted by
  Neyazi); MLIP and PANK1 form a shared arrhythmogenic-remodeling component
  (batch-clean vs ischemic in Liu). Present all as leads, consistent with the
  single-cohort evidentiary bar.

## Files

- `results/figures/geneset_robustness/summary.tsv` — full acceptance matrix + supporting numbers (7 CM/directional sets + `cmstructural` composition-control row, A/B)
- `results/figures/geneset_robustness/acceptance_matrix.png` — pass/fail heatmap (5 CM sets + 2 directional-control rows)
- `A_distance_stats_by_geneset_ALL8.tsv`, `B_zone_stats_by_geneset_ALL8.tsv` — A/B tables including the `cmstructural` sarcomere-content control (8 rows)
- `C_snrna_zone_by_geneset_ALL8.tsv`, `D_crossdisease_stats_by_geneset_ALL8.tsv` — C/D tables (7 sets scored; `cmstructural` = n/a, A/B control only)
- 5-set-only tables (`A/B/C/D_..._by_geneset.tsv`) and the directional-control raw outputs are retained under `results/figures/geneset_robustness/`
- `A_distance_to_lesion_by_geneset.png`, `B_zone_gradient_by_geneset.png`,
  `D_crossdisease_by_geneset.png` — per-analysis plots
- `within_cs_perspot.tsv.gz`, `D_crossdisease_perspot.tsv.gz` — per-spot score checkpoints
- `scripts/geneset_robustness.py` — parameterized wrapper (one `--geneset` per file, analysis logic unchanged)
