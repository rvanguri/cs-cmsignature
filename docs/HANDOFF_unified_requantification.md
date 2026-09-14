# Handoff — Unified re-quantification to remove the cross-cohort batch confound

**Prepared for:** downstream Claude Science agent (HPC/SLURM execution)
**Status:** analysis-motivated engineering spec. No manuscript claim is settled until the acceptance
gate in §7 passes.
**One-line goal:** re-derive the CS-vs-mimic cardiomyocyte signature from a count matrix in which all
cohorts were quantified by *one identical pipeline*, so that technical negative-control genes
(ribosomal/OXPHOS housekeeping) are **not** differentially expressed between CS and mimics — which
they currently are.

---

## 1. Why this is needed (the problem, with evidence)

The committed discovery contrast (`results/de/de_cm/DE_CS_vs_*.tsv`) nominates the four-gene
signature **GJB7, TNNI3K, MLIP, PANK1**. Re-analysis of the *committed pseudobulk table*
(`results/de/cm_pseudobulk_counts.tsv`, 18,790 genes × 142 samples) shows the cross-cohort contrast
is confounded at the **gene level**, not just by design rank-deficiency:

1. **MLIP and PANK1 are not CS-specific vs the arrhythmogenic mimics.** Raw per-sample log2CPM,
   CS vs DCM (Mann-Whitney): MLIP p=0.83 (CS median 11.35 vs DCM 11.26); PANK1 CS *lower* than
   DCM/ARVC. They separate CS only from HCM/NF. Their large "logFC" in the DE table is an artifact of
   the rank-deficient design (the reported `CS_vs_DCM` logFC ≈ raw `CS − HCM`, i.e. the confounded
   contrast collapses onto the lowest-expressing group). The DE table itself flags every gene
   `design_status = rank_deficient_EXPLORATORY_confounded`, `contrast_confounded = True`.

2. **Housekeeping genes score as "CS-specific," even after quantile normalization.** A strict
   pairwise test (CS higher than *every* mimic individually, Δlog2 ≥ 1, one-sided MWU p<0.05) returns
   ~320 genes; **35% are lncRNA/clone loci (RP11-, AC*, LINC*)** and the top protein-coding hits
   include **RPL26, MRPS28, COX16, NDUFC2** (ribosomal / mito-ribosomal / OXPHOS). Forcing every
   sample to an identical distribution by **quantile normalization does not remove them**
   (RPL26 Δ=4.05, MRPS28 Δ=4.18, still "CS-specific"). Ribosomal genes cannot be 16-fold
   disease-specific — this is a per-gene technical offset between cohorts, which QN (a global,
   per-sample transform) cannot fix.

3. **TNNI3K/GJB7 survive the pairwise test but cannot be certified from these data.** TNNI3K is the
   strongest candidate (CS ~64× over every mimic, non-overlapping distributions, protein-coding
   cardiac kinase, effect varies by comparator) but it ranks *below* ribosomal artifacts (11th among
   protein-coding after QN), and the only batch-clean CS contrast available (Liu CS-vs-ICM,
   GSE205734) cannot test it because the 5′ assay barely captures it (baseMean ≈ 50; padj 0.58).

**Root cause.** The unified pipeline (`slurm/02_cellbender → 03_decontx → 05_scanvi → 07_pseudobulk_de`)
harmonizes *ambient signal, droplets, and cell-type labels* — but each cohort was ingested as a
**pre-computed count matrix from its own Cell Ranger run**, with a different reference build and,
critically, a different **intronic-read policy**:

- Neyazi (CS): Cell Ranger 6.1.0, GRCh38, **"including intronic reads"** (per GEO GSM records).
- Reichart (DCM/ARVC/NF): ingested from **cellxgene** (`scripts/download_reichart_cxg.py`) — a
  third-party re-quantification, not the same Cell Ranger run.
- Larson/Chin (HCM/controls): each cohort's own GEO matrix.

snRNA intronic inclusion changes counts most for **long, intron-rich genes and nuclear lncRNAs** —
exactly the classes dominating the artifact list — and shifts the whole library-composition baseline,
which is why even ribosomal genes acquire a spurious CS-vs-mimic offset. Ambient correction and
scANVI cannot undo a quantification-level difference.

> **Bottom line the fix must achieve:** a CS-vs-mimic pseudobulk in which a curated panel of
> housekeeping negative controls shows |log2FC| < 0.5 and no significance. Until that holds, no
> positive CS-specific gene (TNNI3K included) can be trusted from the cross-cohort contrast.

---

## 2. Scope

**In scope:** re-quantify **all snRNA cohorts** (CS + every mimic + controls) from **raw reads**
through **one Cell Ranger configuration** (same version, same reference, `--include-introns=true` for
all), then re-run the existing downstream unchanged. Add a **negative-control QC gate** and re-derive
the signature.

**Out of scope / unchanged:** the Foong Visium arm (`slurm/10_foong_spatial`) and the within-CS
stage/zone analyses — these are within-cohort and **not affected** by the cross-cohort batch. Keep
them as-is. The Liu standalone contrast (`08_liu_standalone`) is already within-study/batch-clean;
keep it as an independent tier but note its capture limitation for TNNI3K/GJB7.

---

## 3. Data sources and access gates (do this first)

| Cohort (code label) | Disease | Accession | Raw-read source | Access |
|---|---|---|---|---|
| Neyazi | CS | GSE319770 / GSE319771 | SRA (SRX… under the GSE) | **Open** |
| Liu | CS, ICM | GSE205734 | SRA | **Open** |
| Larson | HCM | GSE174691 | SRA | **Open** |
| Chin 2022 | HCM, NF | GSE181764 | SRA | **Open** |
| Chin 2021 | NF | GSE161921 | SRA | **Open** |
| Reichart | DCM, ARVC, NF | EGAS00001006374 | **EGA raw FASTQ/CRAM** | **Controlled — EGA DAC / DUA required** |
| (if used) Chaffin/Simonson | DCM/HCM/NF | phs001539 | dbGaP | **Controlled — dbGaP DAR required** |

**Gate A — controlled data.** Reichart (the entire DCM/ARVC arm, and some NF) is EGA-controlled. Raw
re-quantification of that arm **cannot begin** until the EGA Data Access Committee approves a DUA for
`EGAS00001006374`. Start that request immediately; it is the critical-path dependency.

**Gate B — run open arms first.** Neyazi (CS), Liu, Larson, Chin (HCM + NF controls) are SRA-open.
Re-quantify and validate the **open-only** matrix while Gate A is pending; fold Reichart in after
approval **without rework** (the pipeline already supports staged arms, per `docs/pipeline.md` Phase 0).

> If EGA access is denied or too slow, see §8 (fallback): a negative-control-based normalization
> (RUVg/ComBat-seq) on the *existing* matrices, plus reliance on within-study contrasts. It is a
> mitigation, not a cure — the acceptance gate in §7 still governs what may be claimed.

---

## 4. The one canonical quantification config (apply to every sample)

Pin all of the following identically across cohorts. **This uniformity is the entire point** — any
per-cohort deviation reintroduces the confound.

- **Aligner/counter:** Cell Ranger **`cellranger count`**, a single pinned version (recommend the
  current 8.x; whatever is chosen, use it for *all* samples). Record `cellranger --version` in the
  manifest.
- **Reference:** one 10x GRCh38 reference for all (e.g. `refdata-gex-GRCh38-2024-A`). No mixing of
  2020-A / 2024-A. Record the exact reference tarball name + checksum.
- **Intronic reads:** `--include-introns=true` for **every** sample (snRNA standard; also the setting
  Neyazi used). This is the single most important flag — it is the main axis of the current batch.
- **Chemistry:** `--chemistry=auto` unless a sample is known to need an override; log the detected
  chemistry per sample.
- **No per-cohort custom GTF, no exon-only quantification, no third-party re-quant (drop the
  cellxgene ingest for Reichart — go back to EGA FASTQ).**
- Emit the **raw** (unfiltered) feature-barcode matrix per sample so the existing
  `02_cellbender`/`03_decontx` ambient steps run on comparable inputs.

Deliverable of this phase: `$PROJ/raw_requant/<cohort>/<sample>/raw_feature_bc_matrix/` for every
sample, plus `$PROJ/raw_requant/quant_manifest.tsv` (sample, cohort, disease, patient, CR version,
reference, include_introns, chemistry, n_reads, checksum).

### New SLURM step
Add **`slurm/01b_requant.sbatch`** (job array, one `cellranger count` per sample; CPU partition,
high mem, local-SSD scratch for BAM/temp). Insert between `01_download` and `02_cellbender`.
`01_download.sbatch` must be extended to fetch **FASTQ** (SRA via `prefetch`+`fasterq-dump` for open
cohorts; EGA download client for Reichart) rather than processed matrices. Keep the old matrix-ingest
path behind a flag for provenance/comparison, but it must **not** feed the new DE.

---

## 5. Downstream (reuse existing pipeline, no science change)

Once every cohort has raw matrices from the canonical config, run the existing steps **unchanged**:

1. `02_cellbender.sbatch` → per-sample ambient removal (GPU array). Uniform now because inputs are
   uniform.
2. `03_decontx.sbatch` → DecontX integer counts.
3. `04_qc.sbatch` → per-dataset QC h5ad.
4. `05_scanvi.sbatch` (+ `05b_relabel`) → integration + cell-type labels. `study` remains the batch
   key for scANVI (integration is still legitimate; we are removing the *quantification* batch, not
   pretending cohorts are one batch).
5. `06_validate_gates.sbatch` → existing decontam + procurement gates.
6. `pseudobulk_cm.py` (via `07_pseudobulk_de.sbatch`) → new `cm_pseudobulk_counts.tsv` + `_meta.tsv`
   **written to a new dir** (`results/de_requant/…`) so the current committed tables are preserved for
   before/after comparison.
7. `run_pseudobulk_de.R` **and** `run_pseudobulk_de_patientblock.R` → DE on the re-quantified
   pseudobulk. Patient-blocked variant is preferred (Neyazi CS patients contribute multiple
   region/replicate columns).

**Design note (still confounded — be honest):** CS is still a single cohort, so `~0 + disease + study`
stays rank-deficient and `CS_vs_mimic` remains formally confounded with `study`. Re-quantification does
**not** make the design full-rank; it removes the *quantification* component of the batch so that
what's left is closer to biology. The **acceptance gate (§7), not the p-values, decides trust.** Add
the RUVk factors from §6 to the design to absorb residual unwanted variation.

---

## 6. New: negative-control-based residual correction (RUVg)

Even after uniform quantification, add an explicit unwanted-variation control keyed to negative
controls:

- Curate a **technical negative-control gene set** `metadata/neg_control_housekeeping.tsv`: ribosomal
  proteins (`^RPL`, `^RPS`), mito-ribosomal (`^MRPL`, `^MRPS`), core OXPHOS (`NDUF*`, `COX*`, `UQCR*`,
  `ATP5*`), and canonical housekeepers (`ACTB`, `GAPDH`, `B2M`, `TBP`, `PGK1`). These are assumed
  **not** disease-regulated between CS and cardiomyopathy mimics.
- Estimate unwanted-variation factors with **RUVg** (`RUVSeq`) on the pseudobulk using that
  negative-control set; include `k` (start k=1–4, chosen so the negative-control panel flattens —
  see §7) factors as covariates in the limma design.
- Sanity check with **ComBat-seq** (`sva`) as an alternative (batch = `study`) and confirm the two
  give concordant CS-specific gene lists.

Ship this as `scripts/ruv_pseudobulk.R` (input: pseudobulk counts + meta + neg-control list; output:
corrected DE tables + the negative-control diagnostic).

---

## 7. Acceptance gate (this is the deliverable that matters)

Re-quantification + correction is **accepted only if all hold** on the re-derived CS-vs-mimic DE:

**G1 — Negative controls flat (primary gate).** For the housekeeping panel (§6), median |log2FC|
CS-vs-each-mimic < 0.5 and **no** gene FDR<0.05. Report the panel's fold-change distribution
before/after. *Current status: FAILS badly (RPL26/MRPS28/COX16 at ~4 log2). This is the number that
must move.*

**G2 — lncRNA/clone fraction collapses.** Among genes passing the strict pairwise CS-specific test,
the fraction of `RP11-/AC*/LINC*` loci drops from the current 35% toward the transcriptome background
(<~10%).

**G3 — Positive controls behave.** Known CS/granuloma biology (e.g. immune/inflammation markers used
elsewhere in the paper) remains detectable in the appropriate lineages; the CM lineage is not wiped
out.

**G4 — Re-adjudicate the four genes** on the accepted matrix, reporting for each the raw per-mimic
Δlog2 + MWU **and** the corrected DE:
- MLIP, PANK1: expected to remain **non-specific vs DCM/ARVC** (equal). If so, reclassify as a shared
  arrhythmogenic-remodeling program; do not present as CS-specific.
- TNNI3K, GJB7: report whether they survive **once the negative controls are flat**. Only then is
  "CS-specific" defensible for them.

**G5 — Fresh derivation.** Recompute the CS-specific CM gene set from the accepted matrix (pairwise
"higher than every mimic" + RUV-corrected FDR). Deliver the new list; the signature the paper uses
should be defined from *this*, not the confounded table.

Emit `results/de_requant/negative_control_gate.tsv` and a one-page `docs/requant_acceptance.md`
summarizing G1–G5 with before/after numbers.

---

## 8. Fallback if controlled raw data is blocked

If EGA/dbGaP access for Reichart is not obtainable:

- **Do not** re-quantify only the open cohorts and compare them to the still-old Reichart matrix —
  that leaves the DCM/ARVC arm on a different pipeline and reintroduces the confound for the mimics
  that matter most.
- Instead run §6 (RUVg + ComBat-seq negative-control correction) on the **current** matrices as a
  mitigation, and gate on §7-G1 exactly the same way. If G1 cannot be made to pass, the honest
  conclusion is that the cross-cohort four-gene specificity claim is **not supportable** from
  available data, and the paper should stand on the batch-clean tiers (within-CS spatial/stage
  biology + Liu within-study), presenting TNNI3K as a hypothesis-generating candidate for an
  independent CS cohort. (A parallel honest-reframing draft is being prepared separately.)

---

## 9. Concrete task list for the agent

1. [Gate A] File the EGA DAC/DUA request for `EGAS00001006374`; file dbGaP DAR for `phs001539` if
   that cohort will be used. Record status in `docs/data_availability.md`.
2. Extend `scripts/download_geo.py` + `slurm/01_download.sbatch` to fetch **FASTQ** (SRA `prefetch` +
   `fasterq-dump`) for the open cohorts; add an EGA-download path for Reichart (post-approval).
3. Add `slurm/01b_requant.sbatch` running `cellranger count` with the pinned config in §4 (job array,
   one sample per task). Write `raw_requant/quant_manifest.tsv`.
4. Re-run `02→06` unchanged on the re-quantified inputs (new output dirs; do not overwrite committed
   results).
5. Run `pseudobulk_cm.py` → `results/de_requant/cm_pseudobulk_counts.tsv` (+ meta).
6. Add `scripts/ruv_pseudobulk.R` (RUVg + ComBat-seq, §6) and `metadata/neg_control_housekeeping.tsv`.
7. Run `run_pseudobulk_de_patientblock.R` (+ RUV factors) → re-derived DE.
8. Evaluate the **acceptance gate §7**; write `results/de_requant/negative_control_gate.tsv` and
   `docs/requant_acceptance.md`.
9. Re-derive the CS-specific CM signature (G5); hand the new gene list + gate report back for
   manuscript decisions.

## 10. Reproducibility / provenance

- Pin and record: Cell Ranger version, reference tarball + checksum, `--include-introns` value,
  chemistry per sample, and per-cohort read counts in `quant_manifest.tsv`.
- Keep the old (matrix-ingest) results untouched for a documented before/after; the negative-control
  gate §7-G1 is the headline before/after metric.
- All new outputs under `results/de_requant/` and `docs/requant_acceptance.md`; nothing overwrites the
  currently committed tables.

---

### Appendix — key numbers to reproduce (from the committed pseudobulk, pre-fix)

| Gene | CS vs DCM (raw MWU) | CS median log2CPM | DCM | ARVC | HCM | NF | verdict |
|---|---|---|---|---|---|---|---|
| TNNI3K | p=1.7e-16 | 9.99 | 3.45 | 3.92 | 5.78 | 3.75 | candidate CS-specific (unconfirmed) |
| GJB7 | (CS-highest, low abundance) | 3.47 | 2.46 | 2.39 | 0.81 | 1.70 | modest candidate |
| MLIP | p=0.83 | 11.35 | 11.26 | 11.69 | 8.03 | 11.05 | **shared with DCM/ARVC — not CS-specific** |
| PANK1 | CS lower | 6.18 | 6.31 | 6.70 | 4.27 | 5.63 | **shared with DCM/ARVC — not CS-specific** |

Negative-control panel, CS-vs-mimic min Δlog2 **after quantile normalization** (should be ~0; is not):
RPL26 4.05 · MRPS28 4.18 · COX16 3.68 · NDUFC2 3.14 → **primary evidence of gene-level cross-cohort
batch; §7-G1 is the fix's success metric.**

Batch-clean cross-check (Liu CS-vs-ICM, GSE205734, within-study): MLIP padj 0.045, PANK1 padj 0.019
(CS>ICM); TNNI3K padj 0.58, GJB7 padj 0.99 (capture-limited, not evidence of absence).
