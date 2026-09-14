# Handoff — Gene-set robustness of the spatial conclusions

**Prepared for:** downstream Claude Science agent (HPC/SLURM, raw Foong Visium required)
**Standalone.** Does not depend on the re-quantification handoff. This is a self-contained
robustness experiment.
**One-line goal:** determine whether the paper's two spatial conclusions — (i) **cell-autonomy /
anti-paracrine** localization and (ii) **remodeling-stage-independence** — hold when the cardiomyocyte
program is scored on **defensible gene sets** instead of the hand-curated four genes
(GJB7, TNNI3K, MLIP, PANK1).

---

## 1. Why this is needed

The four-gene panel is not objectively justified. Ranked by the pipeline's own `evidence_score`
(`results/de/de_cm/CM_foong_validation.tsv`), the four genes sit at ranks 336 / 420 / 1,978 / 3,781
of ~18,700; the top of that list is entirely different genes (CTNNA3 — a canonical ARVC gene — is
#1 and was not chosen). The panel was selected by a subjective "relevance to arrhythmia" filter, not
by the data. Worse, a pre-specified objective rule shows the cross-cohort "CS-specific" set and the
batch-clean set are **disjoint** (0 of 53 objective genes replicate in the batch-clean Liu CS-vs-ICM
contrast; the only Liu-significant genes, MLIP/PANK1, are the ones that fail specificity).

**Consequence:** the paper's spatial results (cell-autonomy, stage-independence) were all computed on
the four-gene module score. If those conclusions are an artifact of the specific gene choice, the
paper has no spine. If they are **robust across gene sets**, the compartment-level biology stands
independent of which genes are picked, and that becomes the defensible core of the paper. **This
experiment decides that.**

---

## 2. Gene sets to test (provided)

Written to `metadata/` (one symbol per line):

| File | Set | Rationale |
|---|---|---|
| `geneset_original4.txt` | GJB7, TNNI3K, MLIP, PANK1 | the curated panel (baseline / current paper) |
| `geneset_objective53.txt` | 53 genes | pre-specified objective rule: CS > **every** mimic (Δlog2 ≥ 1, one-sided MWU p<0.05) **and** not a technical class (lncRNA/clone, ribosomal, OXPHOS, housekeeping) **and** Foong-replicated |
| `geneset_tnni3k.txt` | TNNI3K only | the single standout gene (largest, most consistent CS-vs-mimic effect) |
| `geneset_mlip_pank1.txt` | MLIP, PANK1 | the only batch-clean (Liu CS>ICM) genes; the shared arrhythmogenic-remodeling pair |
| `geneset_negctrl53.txt` | 53 random CM-expressed genes | **negative control** — must NOT reproduce the spatial pattern (seed=7) |
| `geneset_highconf111.txt` | 111 genes | the pipeline's **batch-clean-concordant** signature (`results/de/tables/high_confidence_cs_set.tsv`): significant + same-sign across every primary contrast **and** Liu. This is the most principled existing signature — but it is **inflammation/fibrosis-dominated** (CD163, CD74, CSF2RA, COL1A1/3A1, FN1…) and contains only PANK1 of the four CM genes. Included as a **contrast set** (see §6). |
| `geneset_inflammation4.txt` | NLRC4, IL1RAP, BACH2, IL7 | the curated granuloma-inflammation axis (paper's positive control). Included as a **directional positive control**: it should localize to **granulomatous** zones, i.e., the *opposite* of the CM candidate sets. |

The objective-53 derivation is fully reproducible from `results/de/cm_pseudobulk_counts.tsv` +
`cm_pseudobulk_meta.tsv` + `CM_foong_validation.tsv`; the rule is stated above and the list is frozen
in the file so results are deterministic.

**Context — the larger signatures already exist.** Claude Science previously produced a 1,196-gene CM
DE pool (`results/de/cs_cm_signature.txt`, the set the four were curated from) and the 111-gene
batch-clean set above. Neither is the paper's headline because the batch-clean set lands on the
already-known granuloma/fibrosis biology, and the four-gene "cardiomyocyte-intrinsic" panel was carved
out of the confounded residual after removing lineage genes. This robustness check is what decides
whether the *spatial compartment-level conclusions* depend on that carve-out.

---

## 3. Data inputs (raw — this is why it needs the cluster)

The committed per-spot table (`results/figures/foong_regional/foong_regional_perspot.tsv.gz`) holds
**only the 4-gene composite score**, so it cannot answer this question. Re-scoring requires the raw
spatial objects:

- **Foong CS Visium** — GEO **GSE314910** (12 sections / 9 patients). Downloader:
  `scripts/download_foong_spatial.sh` (SLURM step 10), `FOONG_GSE=GSE314910`.
- **Kuppe normal controls** — the 4 control sections already used as the scoring baseline.
- **The 4 newly generated arrhythmogenic mimics** — `DCM_LMNA.h5ad`, `LMNA.h5ad`, `ARVC_1.h5ad`,
  `ARVC_2.h5ad` (from `HANDOFF_new_visium_comparators.md`; artifact IDs listed there).

All are **open-access** (no controlled-data gate), so this experiment can run immediately and does not
depend on the EGA/Reichart access that blocks the re-quantification handoff.

---

## 4. Method — re-score, then re-run the identical analyses per gene set

Parameterize the existing scoring so the gene set is an argument; **change nothing else** in the
analysis logic.

1. **Score each gene set** on CM-dominant spots using the same module-score routine already in
   `scripts/foong_regional.py` / `scripts/cm_spatial_crossdisease.py` (Scanpy `score_genes` on
   `argmax lineage == Cardiomyocyte` spots). Add a `--geneset <file>` / `--label <name>` argument;
   loop over the five files in §2. Single-gene sets (TNNI3K) use the gene's log-normalized expression
   directly (score_genes is undefined for n=1).
2. For each gene set, re-run the **three within-CS analyses** and the **cross-disease pattern**,
   exactly as in the published pipeline:

   **(A) Cell-autonomy / anti-paracrine — distance-to-lesion mixed model.**
   Reproduce `foong_distance_stats.tsv`: `mixedlm` of score on within-sample-SD distance-to-lesion,
   spots nested in patients. Report `beta_per_SD`, 95% CI, p, and **the number of patients with
   positive slope (of 8)**. *Published 4-gene result: β=+0.044, CI +0.041/+0.046, 8/8 positive.*
   The **sign and 8/8 consistency** are the claim, not the magnitude.

   **(B) Stage-independence — zone gradient.**
   Reproduce `foong_zone_stats.tsv`: score by zone (preserved / granulomatous / fibrotic),
   patient-level. Report per-patient medians and the paired preserved-vs-fibrotic and
   preserved-vs-granulomatous Wilcoxon. *Published 4-gene result: preserved highest, p≈0.001.*

   **(C) Stage-independence — snRNA cross-check (optional, no raw data needed).**
   On `cm_pseudobulk_counts.tsv` (CS only), score each set and test flatness across Neyazi histologic
   zones (`metadata/neyazi_cs_zones_by_section.tsv`) and vs the Table S1 stage variables if available.
   *Published: flat (Kruskal p=0.35; EF/procedure n.s. after de-confounding).*

   **(D) Cross-disease pattern.**
   Per-gene-set mean score in CM-dominant spots for CS vs Kuppe-normal, vs HCM (Foong-published
   CS-vs-HCM), and vs the 4 arrhythmogenic mimics. *Expected (from 4-gene): CS ≈ ARVC ≈ LMNA > HCM ≈
   normal.* Test whether this ordering is gene-set-dependent.

3. Emit one tidy table per analysis with a `geneset` column, plus a summary matrix (§6).

Suggested new script: `scripts/geneset_robustness.py` (wraps the existing scoring + the three tests);
new SLURM runner `slurm/12_geneset_robustness.sbatch`. Write outputs to
`results/figures/geneset_robustness/` — do not overwrite committed tables.

---

## 5. Reference numbers to reproduce (4-gene baseline)

| Analysis | 4-gene published value |
|---|---|
| (A) distance β per SD | +0.044 (CI +0.041/+0.046), 8/8 patients positive slope, p<1e-260 |
| (B) zone (patient-level) | preserved median ≈ 0.00 > fibrotic ≈ −0.29, granulomatous ≈ −0.30; Wilcoxon p≈0.001 |
| (C) snRNA zone flatness | Kruskal p=0.35; fibrotic-burden Spearman rho=0.17 p=0.46; TX vs LVAD p=0.97 |
| (D) cross-disease | CS ≈ ARVC ≈ LMNA > HCM ≈ normal |

The robustness question is whether (A)–(C) **reproduce in direction and significance** for the
objective-53 and TNNI3K-only sets, and (D)'s ordering is preserved.

**Directional controls (must behave oppositely):** `geneset_inflammation4.txt` and
`geneset_highconf111.txt` are inflammation/fibrosis-weighted, so they should score **highest in
granulomatous/fibrotic zones and near lesions** — the *reverse* of the CM candidate sets. If they do,
the scoring and zone/distance logic are working and the CM sets' preserved/anti-paracrine pattern is
real signal. If the inflammation sets instead look preserved-enriched too, the pipeline is broken and
nothing is interpretable.

---

## 6. Acceptance / interpretation matrix (the deliverable)

Build `results/figures/geneset_robustness/summary.tsv` — rows = gene sets, columns = the four
analyses — and read it as:

- **Robust (paper has a spine):** (A) anti-paracrine sign + ≥7/8 positive slopes, and (B)/(C)
  stage-independence, hold for **original-4 AND objective-53 AND TNNI3K-only**, and the negative
  control does **not** show them. → The cardiomyocyte-intrinsic, stage-independent CS state is a real
  property of the compartment, independent of gene choice. Lead the paper on this; treat specific
  genes as illustrative candidates.
- **Fragile (paper must change):** the conclusions hold for original-4 but **fail** for objective-53
  or TNNI3K-only, or the **negative control also passes**. → The spatial results are an artifact of
  the curated panel (or trivially true for any CM gene set). The signature framing cannot be salvaged;
  the honest paper collapses to "MLIP/PANK1 mark a shared arrhythmogenic-remodeling program
  (batch-clean vs ICM); a broader CS-specific cardiomyocyte signature is not supported."
- **Partial:** report exactly which conclusions survive which sets; the manuscript claims must be
  scoped to the gene sets that pass.
- **Directional controls:** confirm `inflammation4` and `highconf111` score highest in
  granulomatous/fibrotic zones / near lesions (opposite of the CM sets). A pass here validates the
  whole experiment; a failure invalidates all other rows.

A **negative-control pass is as important as the candidate results** — if random CM genes reproduce
the anti-paracrine/preserved-enriched pattern, the pattern reflects CM-content/tissue architecture
(preserved spots are simply cardiomyocyte-richer), not the signature, and (A)/(B) do not support a
program-specific claim at all. Test this explicitly.

---

## 7. Caveats to carry

- Single-gene (TNNI3K) scoring is expression, not `score_genes`; note the methodological difference.
- TNNI3K is a 309 kb gene under-detected by poly-A Visium, so a null spatial result for TNNI3K-only
  may be detection-limited rather than biological — interpret a TNNI3K-only failure with that in mind
  (it does not by itself sink the snRNA finding).
- Mimic N is 1–2 patients; (D) is directional, not powered.
- This experiment tests robustness of the **spatial** conclusions only; it does not resolve the
  cross-cohort specificity confound (that is the separate, and now infeasible, re-quantification
  path).

---

## 8. Task list

1. Add `--geneset/--label` to the scoring path in `scripts/foong_regional.py` (and
   `cm_spatial_crossdisease.py` for §4D); no other logic changes.
2. `scripts/geneset_robustness.py` loops the five `metadata/geneset_*.txt` files over analyses
   (A)–(D); `slurm/12_geneset_robustness.sbatch` runs it against GSE314910 + Kuppe + the 4 mimics.
3. Emit per-analysis tidy tables + `results/figures/geneset_robustness/summary.tsv` and a one-page
   `docs/geneset_robustness_result.md` filling in §6.
4. Report back the summary matrix; manuscript framing decision follows from it.
