# Data and code availability audit — ESC HF Short Report

> Audit of an earlier working tree and draft, kept as written. Since then the Kuppe normal-donor sections have been removed from the analysis, the pipeline and its outputs are committed, and the comparator H&E image for Figure 1A is in `data/eschf/img/`. Current provenance is in `docs/data_availability.md`.

Read-only audit of the `cs-cmsignature` working tree (local HEAD `e376985`, with uncommitted changes) against
`CS_cardiomyocyte_ESCHF_ShortReport_v1.docx`. No analysis was rerun and no code or table was changed. This file is
the only one added. Public record headers were read from GEO (series/sample records, E-utilities), CELLxGENE
(curation API) and GitHub (repository API), as the handoff permits. "Not determinable" means the repository does not record
the fact. The pipeline's raw inputs sit under `/gpfs/data/vangurilab/heartmap-cs/raw` (cluster), which was not accessed.

Line citations are to the current working tree unless prefixed `pre-slim-2026-09-14:`, the git tag holding the retired fetch/QC
pipeline (removed from the tree in `25d9df1`).

## 1. Dataset table

| Dataset | Accession | Source study | Groups and patients used | Used in which reported result(s) | Evidence (script:line) | Status |
|---|---|---|---|---|---|---|
| CS snRNA-seq | GEO GSE319770 | Neyazi et al., *Circulation* 2026 | CS 21 pts / 41 spec (21 of 22 published patients; SNU-1 absent) | all snRNA contrasts (CS vs DCM, ARVC, NF; pooled) → direction agreement, Fig 1C–D, A5b procurement | `pre-slim-2026-09-14:scripts/make_sample_meta.py` (neyazi→GSE319770); `pre-slim-2026-09-14:scripts/build_manifest.py` (raw `*_raw_feature_bc_matrix.h5` → CellBender); `scripts/eschf/00_inputs.py` (pseudobulk from cluster dir); `scripts/eschf/a5b_procurement.R:16–19,47` | confirmed |
| CS snRNA clinical table | Neyazi et al. Supplementary Table S1 (sheet "Sarcoid") | Neyazi et al. 2026 | 21 CS patients (Tissue Source; LVEF) | A5b procurement (VAD vs explant), clinical columns | `scripts/eschf/00_inputs.py` (Table S1 read → `data/external/neyazi_clinical.tsv`); `a5b_procurement.R:20–22` | confirmed |
| CS Xenium | GEO GSE319771 | Neyazi et al. 2026 | none | none | no reference in any script in tree or at `pre-slim-2026-09-14` (only in docs) | confirmed **not read** |
| DCM/ARVC/NF snRNA-seq | CELLxGENE collection e75342a8-0f3b-4ec5-8ee1-245a23e0f7cb, dataset 65badd7a-9262-4fd1-9ce2-eb5dc0ca8039 (version 6d0b8b8a-2b22-4578-b8ee-8ca4ab893594, "All cells"); raw data EGA EGAS00001006374 | Reichart et al., *Science* 2022 | ARVC 8 pts / 8 spec; DCM 52 pts / 52 spec; NF 18 pts / 18 spec | all snRNA contrasts; the CS-vs-NF contrast (18 of 24 NF donors); cell-type label reference for scANVI | `pre-slim-2026-09-14:docs/pipeline.md` (CELLxGENE "All cells" URL, raw counts in `raw.X`); `pre-slim-2026-09-14:scripts/download_reichart_cxg.py`; `pre-slim-2026-09-14:scripts/run_scanvi.py` (Reichart labels as reference) | confirmed from code/docs; the exact file on the cluster was not verified. **EGA itself was never accessed** (processed CELLxGENE release only) |
| HCM + NF snRNA-seq | GEO GSE174691 | Larson/Codden/Chin (Chin lab), 2021 | HCM 9 pts / 9 spec | **in the model fit** (shares variance/filtering); HCM contrast not reported | patient IDs = GSM5322568–76, all 9 in GSE174691; `a5b_procurement.R:25,31,47` fits all cohort patients; contrasts at `:36` are DCM, ARVC and NF only | confirmed |
| HCM + NF snRNA-seq | GEO GSE181764 | Codden/Chin, 2022 | HCM 7 pts / 7 spec; NF 2 pts / 2 spec | in the model fit; the 2 NF donors enter CS vs NF | patient IDs = GSM5511023–31, all 9 in GSE181764 | confirmed |
| NF snRNA-seq | GEO GSE161921 | Larson/Chin, 2021 (IVS) | NF 4 pts / 4 spec | CS vs NF (4 of 24 NF donors) | patient IDs = GSM4923401–04, all 4 in GSE161921 | confirmed |
| 5′ CS/control snRNA | GEO GSE205734 | Liu et al. 2022 | none in the analysed cohort | none | `pre-slim-2026-09-14:scripts/make_sample_meta.py` and `pre-slim-2026-09-14:slurm/04_qc.sbatch` processed it; `run_scanvi.py` keeps it standalone; absent from `snrna_patient_cohort.tsv` | confirmed **not used** |
| CS Visium | GEO GSE314910 | Foong et al., *J Card Fail Intersect* 2026 | 11 of 12 deposited sections (CS_109 excluded), 8 of 9 patients | all spatial results (Fig 1A–B, spatial DE, distance model, pathway scores) | `scripts/spatial_first/00_export_from_raw.py:15,52–56`; `scripts/eschf/a2a6_spot_export.py:43,85–95`; `scripts/eschf/a0_primary_export.py` (load_cs); `data/spatial_first/spatial_first_config.json` (11 analysed) | confirmed; cluster copy vs GEO identity **not determinable** (§2.1) |
| Comparator Visium | none (pending) | this study | 6 sections, 6 patients: ARVC ×2, LMNA DCM ×2, HCM, Chagas | all spatial comparisons | `00_export_from_raw.py` (`COMP_DIR`); local copy `~/Downloads/new visium` | confirmed in use; **not deposited** |
| Non-failing Visium (normal) | CELLxGENE collection 8191c283-0816-424b-9b61-c3e1d6258a77, 4 "normal" Visium datasets (§2.5) | Kuppe et al., *Nature* 2022 | 4 sections (NORMAL_5262ff55, _661054b2, _6f7b96ec, _889e6502) | enter the **pooled per-gene z-transform** that underlies argmax lineage assignment for every spot (CS and comparator), i.e. every spatial result. No reported contrast uses them as a group | `00_export_from_raw.py` (load_normal, `kuppe_control/*.h5ad`; pooled z over all sections before argmax); `a2a6_spot_export.py` and `a0_primary_export.py` (same); `data/spatial_first/meta_panel.tsv` | confirmed (file prefixes match CELLxGENE dataset version IDs) |
| Published cardiomyocyte DE | Chaffin et al. *Nature* 2022, Supplementary Tables 6 (DCM vs NF) and 7 (HCM vs NF) | Chaffin et al. 2022 | Cardiomyocyte rows only | Results "lower in disease" counts and backgrounds; Fig 1D columns 4–5 | `scripts/spatial_first/06_published_reference.py:36,39–53` (filter `Cell Type == Cardiomyocyte`; columns Gene, Ensembl ID, `CellBender: logFC`, `CellBender: Adjusted P-Value`, `Background Contamination Flag`); `scripts/eschf/a3a4_qc_nonfailing.py:83,119` | confirmed. Raw data (dbGaP/Broad portal) never accessed |
| Published proteomics | Hamed et al. 2026, ESM 4, 5, 7 | Hamed et al. 2026 | — | **no reported number**. Read by `scripts/eschf/a5_hamed.py` (A5); none of its output values (overlap n, rho, counts) occur in the manuscript text or figure | `a5_hamed.py`; harness A5 block; text search of the draft | confirmed: interpretive citation only |
| KEGG pathway gene sets | hsa00190, hsa00020, hsa00071 (also hsa01100, hsa04512) | KEGG REST `rest.kegg.jp link/hsa` | — | pathway levels (A10/A14), Fig 1B, Results pathway counts | `metadata/geneset_kegg_*.txt` header lines: hsa00190 retrieved 2026-10-05; hsa00020 and hsa00071 retrieved 2026-10-06. KEGG **release number not recorded**. Read by `scripts/eschf/a10_pathway_levels.py` | confirmed (date); release not determinable |
| Gene annotation (locus QC) | Ensembl REST `lookup/symbol/homo_sapiens`; UCSC hg38 `genomicSuperDups` | Ensembl; UCSC | — | qc_pass (gene span, segdup) → QC-passing subsets, Fig 1D marks | `scripts/spatial_first/03_locus_qc.py:52` (Ensembl URL), cached `ensembl_annotation_tested.tsv`, `segdups_hg38.tsv.gz` (`--refresh` refetches) | source confirmed; **Ensembl release and fetch date not recorded** |

Notes on §1 items the handoff asked for explicitly:

- **HCM in the fit.** `a5b_procurement.R:47` calls `fit_contrasts(coh$patient, "full")` on all 121 patients in
  `snrna_patient_cohort.tsv`, so `disease` includes the HCM level. The design is `~0 + disease + anatomy3` (`:31`),
  `filterByExpr(group = disease)` (`:29`), and the contrasts are DCM, ARVC and NF only (`:36`). Lines `:48–52`
  assert that this refit reproduces the committed `snrna_CS_vs_{DCM,ARVC,NF}.tsv` to within 1e-8. The 16 HCM
  patients therefore contribute to gene filtering, TMM normalisation, voom weights and the eBayes prior. The
  manuscript's "HCM was excluded" is accurate for the contrast but not for the fit.
- **Pooled CS vs DCM+ARVC** is inverse-variance pooling of the two committed contrasts (`scripts/spatial_first/04_snrna_pool.py:57`, `pool_snrna`); no separate fit.
- **NF donors (24):** Reichart 18, Chin 2022 2, Chin 2021 4.
- **Reichart / Litviňuková.** The 18 Reichart NF donor IDs are D1, D2, D4–D7, H2–H7, H46, H49, H51, H53, H55, H67. Twelve
  of them (D1, D2, D4–D7, H2–H7) also appear as donor IDs in the Litviňuková et al. Heart Cell Atlas collection
  on CELLxGENE (b52eb423-5d0d-4645-b217-e1c6d38b2e72). Matching ID strings across two collections from the same group indicate
  shared donors. Whether the Reichart nuclei for those donors are the Litviňuková data (re-used) or newly generated
  is **not determinable** from the repository. Whether EGAS00001006374 holds their raw data was not checked (EGA not accessed).
- **CS patient count.** 21 of 22: SNU-1 is listed in Neyazi Table S1, but none of the 41 GSE319770 samples belongs to it
  (`results/eschf/a11_snrna_cohort_accounting.tsv`, `a11_raw_donor_check.tsv`; reconfirmed against the 41 GEO sample titles).
- **Processing route per GEO dataset** (from `pre-slim-2026-09-14:scripts/build_manifest.py`). Neyazi deposits raw
  `*_raw_feature_bc_matrix.h5` → CellBender. Larson/Chin deposits are filtered `matrix.mtx/features/barcodes`
  → DecontX only. Reichart CELLxGENE `raw.X` → DecontX only (`pre-slim-2026-09-14:docs/pipeline.md`). Which exact files reached the
  cluster is not recorded beyond `pre-slim-2026-09-14:scripts/download_geo.py` (whole GEO supplementary directory per series).

## 2. Spatial transcriptomics

### 2.1 CS (GSE314910)

| Section | GSM | Analysed | Probe set (h5 header) | Space Ranger (h5 header) | Barcodes |
|---|---|---|---|---|---|
| CS_100 | GSM9416128 | yes | v2.0 | noversion | 3740 |
| CS_101-1 | GSM9416129 | yes | v2.0 | spaceranger-3.1.2 | 4137 |
| CS_101-2 | GSM9416130 | yes | v2.0 | spaceranger-3.1.2 | 4302 |
| CS_102-1 | GSM9416131 | yes | v2.0 | spaceranger-3.1.2 | 4090 |
| CS_102-2 | GSM9416132 | yes | v2.0 | spaceranger-3.1.2 | 1800 |
| CS_103 | GSM9416133 | yes | v2.0 | spaceranger-2.1.1 | 4974 |
| CS_104 | GSM9416134 | yes | v2.0 | noversion | 3676 |
| CS_106 | GSM9416135 | yes | v2.0 | noversion | 908 |
| CS_107 | GSM9416136 | yes | v2.0 | spaceranger-2.1.1 | 4571 |
| CS_108-1 | GSM9416137 | yes | v2.0 | spaceranger-2.1.1 | 4976 |
| CS_108-2 | GSM9416138 | yes | v2.0 | spaceranger-2.1.1 | 4709 |
| CS_109 | GSM9416139 | no (excluded) | v2.0 | noversion | 864 |

- Twelve sections from nine patients are deposited; eleven sections from eight patients are analysed.
- **The analysed files are a cluster copy.** The loaders glob `{RAW}/foong_spatial/**/*filtered_feature_bc_matrix*.h5` with
  `RAW=/gpfs/data/vangurilab/heartmap-cs/raw` (`00_export_from_raw.py:15,52–56`; `a2a6_spot_export.py:43,85–95`), and parse a
  GSM prefix from each filename. The copy therefore uses GEO's flat `GSM…_` naming. No checksum or file size is recorded, so
  **identity with the GEO deposit is not determinable** without cluster access. The headers above are from the GEO copies.
- **GEO metadata mismatch.** Every GSE314910 sample record states "Space Ranger Software (version 3.1.2)". The h5
  headers say 3.1.2 for four sections, 2.1.1 for four and `noversion` for four. Chemistry (Visium V4 Slide) and probe set
  (v2.0) are uniform. GEO records "Archival FFPE tissue" and CytAssist preparation (CG000518).

### 2.2 Comparators (6 sections) — deposition status

Not deposited. The draft says "[PENDING]", and the repo's `docs/data_availability.md` says "pending GEO deposition".
That edit exists on GitHub only (§4.2). Each comparator folder in `~/Downloads/new visium` holds four files:
`filtered_feature_bc_matrix.h5`, `scalefactors_json.json`, `tissue_hires_image.png` and `tissue_positions.csv`.
None of these is in the repository. A GEO submission would also need the following, none of which is in the repository or
that folder:

1. raw FASTQ per library;
2. raw/unfiltered matrix and the remaining Space Ranger `spatial/` outputs (low-res image, CytAssist image, aligned fiducials), plus
   web summaries / `metrics_summary.csv`;
3. the original H&E (or IF) microscope image;
4. per-sample metadata: diagnosis, age/sex, procurement (explant/VAD/biopsy/autopsy), region, preservation, probe set, slide/area ID;
5. a protocol and data-processing description, with the two Space Ranger versions stated per sample;
6. ethics/consent wording for deposition.

### 2.3 Tissue preservation (comparators)

**Not recorded** anywhere in the repository or the local comparator folder (no Space Ranger config, manifest or
sample sheet). The repo docs say only "probe-based". The probe-based chemistry fits either FFPE or fixed-frozen, so
the preservation cannot be inferred from it.

### 2.4 Probe set and Space Ranger per comparator

| Section | Panel label (`meta_panel.tsv`) | Probe set | Space Ranger | Barcodes |
|---|---|---|---|---|
| ARVC_1 | panelB | v2.1.0 | 4.1.0 | 3579 |
| ARVC_2 | panelB | v2.1.0 | 4.1.0 | 3186 |
| Chagas | panelA | v2.0 | spaceranger-2.1.1 | 4691 |
| DC LMNA | panelB | v2.1.0 | 4.1.0 | 4199 |
| HCM | panelA | v2.0 | spaceranger-2.1.1 | 4165 |
| LMNA | panelB | v2.1.0 | 4.1.0 | 4970 |

The comparators span two probe sets and two Space Ranger versions; the CS sections use v2.0 throughout. No GEO record exists
yet to compare against.

### 2.5 Kuppe normal sections

The section names `NORMAL_<prefix>` in `meta_panel.tsv` match the CELLxGENE `dataset_version_id` prefixes of four Visium
"normal" datasets in collection 8191c283-0816-424b-9b61-c3e1d6258a77:
- NORMAL_5262ff55: dataset 77e04629-a8fa-4d34-bf5f-622d9aed9f35
- NORMAL_889e6502: dataset 58af04fb-de1f-4ba3-8d19-d26610526f14
- NORMAL_661054b2: dataset 4fd2ee79-ab3a-4827-a773-1b7dcd099307
- NORMAL_6f7b96ec: dataset 2257d8b0-7385-48b8-8dca-c5772592c5ea

No download script for these files exists in the tree or at the tag. The match is by version-ID prefix of the filenames.

## 3. Published results used as inputs

1. **Chaffin 2022.** ST6 (`2021-02-03277C-ST6.DCMvsNF.xlsx`) and ST7 (HCM vs NF) are distilled to
   `data/published/chaffin2022_cardiomyocyte_{DCM,HCM}vsNF.tsv` (`06_published_reference.py:36,39–53`; columns above).
   They enter the Results counts of genes lower in disease vs NF, their binomial backgrounds, and Figure 1D columns 4–5
   (`a3a4_qc_nonfailing.py:83,119` → `a4_published_nf.tsv`). The distilled TSVs are tracked in git (redistribution: §5).
2. **Hamed 2026.** MOESM4, MOESM5 and MOESM7 are read by `a5_hamed.py`. **No number now in the manuscript or figure
   depends on them.** The A5 outputs (overlap count, correlation, direction counts) do not appear in the draft text. Hamed is
   cited only for interpretation.
3. **KEGG.** Source `rest.kegg.jp link/hsa`, stored as `metadata/geneset_kegg_*.txt`. Retrieved 2026-10-05 (hsa00190,
   hsa01100, hsa04512) and 2026-10-06 (hsa00020, hsa00071). The KEGG release number is not recorded. `metadata/` is untracked.
4. **Locus annotation.** Ensembl REST symbol lookup (`03_locus_qc.py:52`): the release was whatever was current at fetch and
   is **not recorded**, nor is the fetch date. UCSC `genomicSuperDups` for hg38 is cached as `segdups_hg38.tsv.gz`; its table
   date is not recorded.

## 4. Code

1. **Repository.** https://github.com/rvanguri/cs-cmsignature: public, MIT licence (`LICENSE` present locally and
   detected by GitHub), default branch `main`. The name refers to the withdrawn signature analysis. Proposed rename:
   `cs-cardiomyocyte-transcription` (or `cs-noninflamed-myocardium`). GitHub redirects the old URL after a rename.
   **Not renamed.**
2. **No commit or tag produces the manuscript numbers.**
   - Local HEAD `e376985` ("revised analysis and results", 2026-09-14) is the last local commit.
   - GitHub `main` is **5 commits ahead** (`45d1a8c`), with web edits not in the local clone:
     - README trimmed;
     - `docs/DISCREPANCIES.md` deleted;
     - `docs/data_availability.md` "pending" → "pending GEO deposition";
     - the `scripts/README.md` paragraph pointing to tag `pre-slim-2026-09-14` removed;
     - `00_export_from_raw.py` `RAW` path blanked to `""`.
   - The only tag is `pre-slim-2026-09-14`.
   - The ESC HF pipeline is **untracked**: `scripts/eschf/`, `results/eschf/`, `data/eschf/`, `data/external/`, `metadata/`,
     the figures and `data/snrna/snrna_cm_patient_counts.tsv.gz`.
   - Modified but uncommitted: `scripts/manuscript_results.py`, `results/manuscript_results.tsv`, `DISCREPANCIES.md` and
     `scripts/spatial_first/04_snrna_pool.py`.
   - Pushing requires reconciling with the five GitHub-only commits first. In particular, the GitHub deletion of
     `docs/DISCREPANCIES.md` and the blanked `RAW` path conflict with the local tree.
3. **Harness.** `scripts/manuscript_results.py` is in the tree (modified, uncommitted). Its last written output,
   `results/manuscript_results.tsv` (mtime 2026-10-06 15:29), records 870 OK and 2 NOT_REPRODUCIBLE, 0 mismatched. It was **not rerun** for this audit.
   The figure v5 session reported the same totals.
4. **Reported numbers that depend on files outside the repository:**
   - cluster `RAW/foong_spatial`, `RAW/comparators`, `RAW/kuppe_control` (spatial export; the committed exports in
     `data/eschf/` and `data/spatial_first/` are what the analysis starts from, but `data/eschf/` is untracked);
   - the cluster snRNA pseudobulk directory read by `scripts/eschf/00_inputs.py` (→ `snrna_cm_patient_counts.tsv.gz`, untracked);
   - `data/external/` (Neyazi Table S1/S2 xlsx; Hamed ESM xlsx), untracked;
   - `metadata/` (KEGG sets, structural/ECM gene sets, tested-gene universe), untracked;
   - the comparator Visium folder (local only).

## 5. Draft Data and Code Availability statement (confirmed usage only)

> **Data availability.**
> - Visium data for the cardiac sarcoidosis sections are available from GEO (GSE314910; 11 of the 12 deposited sections
>   analysed)[ref]. Visium data for the six comparator sections will be deposited in GEO [ACCESSION PENDING].
> - Non-failing myocardial Visium sections (four "normal" datasets) were obtained from CELLxGENE (Kuppe et al.)[ref] and were
>   used to standardise lineage scores.
> - Single-nucleus RNA-seq data were obtained from GEO for cardiac sarcoidosis (GSE319770)[ref], hypertrophic cardiomyopathy
>   and non-failing donors (GSE174691, GSE181764)[ref], and non-failing donors (GSE161921)[ref].
> - Dilated and arrhythmogenic cardiomyopathy and non-failing donors were obtained from the CELLxGENE release of Reichart et al.
>   (collection e75342a8-0f3b-4ec5-8ee1-245a23e0f7cb; raw data under EGA accession EGAS00001006374)[ref].
> - Clinical variables for the cardiac sarcoidosis donors were taken from Supplementary Table S1 of Neyazi et al.[2].
>   Cardiomyocyte differential-expression results versus non-failing myocardium were taken from Supplementary Tables 6 and 7
>   of Chaffin et al.[6].
> - KEGG pathway gene sets were retrieved from the KEGG REST service on 5–6 October 2026.
>
> **Code availability.** Analysis code, derived tables and the claim-verification script are available at
> https://github.com/rvanguri/cs-cmsignature (MIT licence) [TAG/COMMIT PENDING].

## 6. Draft [dataset] references

Authors are the first six GEO contributors per record. For the CELLxGENE entries they are the source publication's authors, as
CELLxGENE lists no separate dataset authors. The year is the public release date of the record.

- [dataset] Foong YH, Singh A, Malik M, Zhou F, Liu C, Willet D, et al., 2026, Spatial transcriptomic signatures of cardiac sarcoidosis, Gene Expression Omnibus, GSE314910.
- [dataset] Neyazi M, Venturini G, Brown KJ, Choi Y, Gorham JM, Layton OG, et al., 2026, Integrative Molecular Analyses of Inflammatory and Autoimmune Signals in Cardiac Sarcoidosis [snRNA-Seq], Gene Expression Omnibus, GSE319770.
- [dataset] Larson A, Codden CJ, Chin MT, 2021, Altered Intercellular Communication and Extracellular Matrix Signaling as a Disease Mechanism in Human Hypertrophic Cardiomyopathy, Gene Expression Omnibus, GSE174691.
- [dataset] Codden CJ, Chin MT, 2022, Human Obstructive and Non-obstructive Hypertrophic Cardiomyopathy Disease Mechanisms, Gene Expression Omnibus, GSE181764.
- [dataset] Larson A, Chin MT, 2021, A method for cryopreservation and single nucleus RNA-sequencing of normal adult human interventricular septum heart tissue reveals cellular diversity and function, Gene Expression Omnibus, GSE161921.
- [dataset] Reichart D, Lindberg EL, Maatz H, Miranda AMA, Viveiros A, Shvetsov N, et al., 2022, DCM/ACM heart cell atlas: All cells (Pathogenic variants damage cell composition and single cell transcription in cardiomyopathies), CELLxGENE Discover, collection e75342a8-0f3b-4ec5-8ee1-245a23e0f7cb, dataset 65badd7a-9262-4fd1-9ce2-eb5dc0ca8039.
- [dataset] [AUTHORS — EGA record not read], [YEAR], [TITLE], European Genome-phenome Archive, EGAS00001006374.
- [dataset] Kuppe C, Ramirez Flores RO, Li Z, Hayat S, Levinson RT, Liao X, et al., 2022, Spatial multi-omic map of human myocardial infarction (Visium, normal myocardium), CELLxGENE Discover, collection 8191c283-0816-424b-9b61-c3e1d6258a77, datasets 77e04629-a8fa-4d34-bf5f-622d9aed9f35, 58af04fb-de1f-4ba3-8d19-d26610526f14, 4fd2ee79-ab3a-4827-a773-1b7dcd099307, 2257d8b0-7385-48b8-8dca-c5772592c5ea.

## 7. Datasets listed vs used

**Listed in the draft but not used as listed.**

- EGAS00001006374 is given as the Reichart source, but the pipeline used the CELLxGENE processed release; EGA was never accessed.
- GSE174691 and GSE181764 are both labelled hypertrophic cardiomyopathy, but GSE181764 also supplies 2 of the 24 NF donors.

GSE319771 (Xenium) is not in the draft and was not read.

**Used but not listed in the draft.**

- Kuppe et al. CELLxGENE normal Visium (4 sections; enters lineage standardisation for every spot). Kuppe is not in the reference list.
- The Reichart CELLxGENE collection (the actual source).
- Neyazi Supplementary Table S1 (clinical; covered by reference 2, but not mentioned as a data source).
- Chaffin Supplementary Tables 6 and 7 (reference 6 is cited; not mentioned in the statement).
- KEGG (with retrieval date); Ensembl REST / UCSC genomicSuperDups (locus QC).
- Possibly Litviňuková et al. 2020 (12 shared NF donor IDs; §1).

**Downloaded/processed but not used in any reported number:** GSE205734 (Liu, 5′), GSE319771 and the Hamed ESM files.

## 8. Open items for Rami

1. **Comparator tissue preservation** (FFPE vs fixed-frozen): not recorded anywhere; needed for Methods and GEO.
2. **Comparator deposition:**
   - FASTQ and full Space Ranger outputs;
   - original histology images;
   - per-sample clinical metadata;
   - consent wording for public deposit;
   - whether to deposit both Space Ranger versions as run or re-run under one version.
3. **GSE314910 metadata:** GEO says Space Ranger 3.1.2 for all; headers show 3.1.2/2.1.1/noversion. Request a correction from the depositors, or describe it as is.
4. **HCM wording:** HCM patients are in the fit, though their contrast is not reported. The current "HCM was excluded" sentence needs rewording.
5. **Reichart vs Litviňuková:** confirm from the Reichart paper whether the 12 shared NF donors are re-used Heart Cell Atlas data. If so, whether to cite Litviňuková et al. (*Nature* 2020;588:466–472).
6. **Kuppe:** add the citation and [dataset] entry; decide how to describe the role (lineage-score standardisation).
7. **Repository:**
   - reconcile the local tree with the 5 GitHub-only commits;
   - commit the ESC HF pipeline (`scripts/eschf/`, `metadata/`, `data/eschf/`, derived tables);
   - tag the commit that produces the submitted numbers;
   - decide on the rename;
   - decide whether to restore the blanked `RAW` path or document it as cluster-only.
8. **Redistribution:** whether committing `data/external/` (publisher supplementary xlsx for Neyazi, Hamed) and the tracked Chaffin distillations is permitted. The alternative is fetch scripts.
9. **Versions not recorded:** the KEGG release and the Ensembl release used for locus QC. State the retrieval date only, or re-fetch with version capture (a rerun, not done here).
10. **EGA record** (EGAS00001006374): the title/authors placeholder needs filling from the EGA study page (not read).
