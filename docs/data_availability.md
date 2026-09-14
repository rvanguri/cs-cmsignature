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

- **Spatial arm.** `data/spatial_first/pseudobulk_T2000.tsv.gz` and its companion detection and
  metadata tables are per-section, cardiomyocyte-enriched pseudobulk exported from the raw Visium
  objects (CS from GSE314910, comparators from this study) by
  `scripts/spatial_first/00_export_from_raw.py`. That stage needs the raw objects and is therefore
  **not** part of `run_all.sh`; the pipeline starts from the committed export.
- **Single-nucleus arm.** `data/snrna/snrna_CS_vs_*.tsv` are per-gene contrast summaries computed
  from the snRNA-seq cohorts above (cardiomyocyte nuclei, patient-level pseudobulk);
  `data/snrna/snrna_patient_cohort.tsv` lists the contributing patients by study and disease.
- **Published reference.** `data/published/chaffin2022_cardiomyocyte_*.tsv` are the cardiomyocyte
  DCM-vs-NF and HCM-vs-NF tables transcribed from the Chaffin supplement.
  `scripts/spatial_first/06_published_reference.py --src-dir <dir>` re-derives them from the
  original Excel files if needed.
- Fetch and preprocessing code for the raw single-nucleus matrices is not in this tree; it is in
  git history at tag `pre-slim-2026-09-14`.

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
