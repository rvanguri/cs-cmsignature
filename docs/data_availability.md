# Data availability

All datasets used in the manuscript are de-identified. No raw data is redistributed in this
repository; it carries analysis code, the derived inputs the pipeline starts from
(`data/`), and the result tables and figure it produces (`results/`, `figures/`).

## Datasets used by this repository

| Cohort (label in code) | Modality | Disease(s) | Repository | Accession | Source publication |
|---|---|---|---|---|---|
| This study (comparators) | spatial transcriptomics (Visium, probe-based) | ARVC, LMNA dilated, hypertrophic and Chagas cardiomyopathy | — | pending GEO deposition | this manuscript |
| Foong | spatial transcriptomics (Visium, probe-based) | Cardiac sarcoidosis | GEO | GSE314910 | Foong et al., *J Card Fail Intersect* 2026;2:278-283 |
| Neyazi | single-nucleus RNA-seq | Cardiac sarcoidosis | GEO | GSE319770 / GSE319771 | Neyazi et al., *Circulation* 2026;153:2011-2028 |
| Reichart | single-nucleus RNA-seq | Dilated & arrhythmogenic cardiomyopathy; non-failing control | EGA (open distribution via CELLxGENE) | EGAS00001006374 | Reichart et al., *Science* 2022;377:eabo1984 |
| Larson | single-nucleus RNA-seq | Hypertrophic cardiomyopathy | GEO | GSE174691 | Larson et al., *Sci Rep* 2022, [10.1038/s41598-022-08561-x](https://doi.org/10.1038/s41598-022-08561-x) |
| Chin (2022) | single-nucleus RNA-seq | Hypertrophic cardiomyopathy (+ non-failing) | GEO | GSE181764 | Codden et al., *Int J Mol Sci* 2022, [10.3390/ijms23020946](https://doi.org/10.3390/ijms23020946) |
| Chin (2021) | single-nucleus RNA-seq | Non-failing control | GEO | GSE161921 | Larson et al., *BMC Med Genomics* 2021, [10.1186/s12920-021-01011-z](https://doi.org/10.1186/s12920-021-01011-z) |
| Chaffin | published cardiomyocyte differential expression tables | DCM vs non-failing; HCM vs non-failing | publication supplement | — | Chaffin et al., *Nature* 2022;608:174-180 |

## What is committed, and where it came from

- **Spatial exports.** Exported on the cluster from the raw Visium objects (CS from GSE314910, comparators from this
  study) by the three export stages. These stages need the raw objects and are not part of
  `run_pipeline.sh all`, which starts from the committed exports.
  - `scripts/spatial_first/00_export_from_raw.py` writes `data/spatial_first/`: per-sample cardiomyocyte pseudobulk
    beyond 0.5 mm of immune-dominant spots, detection fractions, section metadata, region selection.
  - `scripts/eschf/a0_primary_export.py` writes `data/eschf/primary_v2/`: the same over all cardiomyocyte-dominant
    spots, with spot counts. `a0_export_log.txt` is that run's log.
  - `scripts/eschf/a2a6_spot_export.py` writes `data/eschf/primary_v2/eschf_perspot_T2000.tsv.gz`: per-spot lineage,
    program scores and distance to the nearest immune-dominant spot.
  No normal-donor sections are included.
- **H&E images.** `data/eschf/img/*.jpg` are the Visium high-resolution tissue images, downscaled by 0.6 (JPEG
  quality 90). Their scale factors are in `image_meta.json`. The CS images come from GSE314910. `ARVC_1.jpg` is the
  comparator section shown in Figure 1A, from this study's data.
- **Annotation.** `data/spatial_first/ensembl_annotation_tested.tsv`, `gene_span.tsv` and `segdups_hg38.tsv.gz`
  are cached Ensembl and UCSC lookups. `results/spatial_first/q2_locus_qc_ensembl_fetched.tsv` and
  `q2_locus_qc_ensembl_unresolved.txt` cover the remaining tested genes, so locus QC runs offline.
- **Single-nucleus contrasts.** `data/snrna/snrna_CS_vs_*.tsv` are per-gene contrast summaries computed from the
  snRNA-seq cohorts above (cardiomyocyte nuclei, patient-level pseudobulk). `snrna_cm_patient_counts.tsv.gz`
  holds the patient pseudobulk counts the procurement and per-study analyses refit, and
  `snrna_patient_cohort.tsv` lists the contributing patients by study and disease.
- **Publisher supplements (not redistributed; download before running).** `run_pipeline.sh` checks for these files and stops
  with this list if any is missing:
  - `data/external/neyazi_table_s1.xlsx` and `neyazi_table_s2.xlsx`: Supplemental Tables S1 and S2 of Neyazi et al.,
    *Circulation* 2026, doi:10.1161/CIRCULATIONAHA.126.079304, saved under these names. S1 must keep its `Sarcoid` sheet.
    They are read by `00_inputs.py` (which writes `neyazi_clinical.tsv`), `a11_*` and the harness.
  - `data/external/44161_2026_880_MOESM4_ESM.xlsx`, `MOESM5`, `MOESM7`: supplementary files of Hamed et al.,
    *Nat Cardiovasc Res* 2026, doi:10.1038/s44161-026-00880-w, under the publisher's download names. They are read
    by `a5_hamed.py`. No number in the short report depends on them.
  `data/external/neyazi_geo_snrna_sample_listing.txt` (the GEO sample listing for the Neyazi snRNA series) is committed.
- **Published reference.** `data/published/chaffin2022_cardiomyocyte_*.tsv` are the cardiomyocyte DCM-vs-NF and
  HCM-vs-NF tables transcribed from the Chaffin supplement. `scripts/spatial_first/06_published_reference.py
  --src-dir <dir>` re-derives them from the original Excel files.
- Fetch and preprocessing code for the raw single-nucleus matrices is not in this tree. It is in git history at tag
  `pre-slim-2026-09-14`.

## Constraints on use

- **The spatial and single-nucleus arms share neither CS cases nor comparators.** The
  single-nucleus arm is a direction check on independent data, not a replication of effect size;
  the two arms are not pooled.
- **Cross-cohort single-nucleus contrasts are confounded.** CS and each comparator come from
  different studies, so study and disease are not separable in those contrasts. They are read for
  direction only, and `results/spatial_first/snrna_gene_length_regression.tsv` documents the
  gene-length trend that excludes CS-vs-HCM from the pooled comparator.
- **Stage heterogeneity.** The single-nucleus CS cohort spans VAD, transplant and autopsy/myectomy
  procurement, whereas all Visium sections (CS and comparators) are end-stage transplant.
- **Platform.** Visium here is probe-based; the single-nucleus and published reference data are
  polyA. Effect sizes are not comparable across those platforms, which is why the non-failing
  comparison is drawn from published tables rather than recomputed.
- **De-identification.** All public samples were de-identified by the original data generators
  before deposition. Consult each accession's record and publication for consent and data-use terms
  before reuse.
