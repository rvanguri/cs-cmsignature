# Data availability

All datasets used in this study are publicly available and de-identified. No raw data is
redistributed in this repository; only analysis code and derived result tables are included. Each
dataset is fetched from its public accession by the step-01 download scripts
(`scripts/download_geo.py`, `scripts/download_reichart_cxg.py`) and, for the spatial series,
`scripts/download_foong_spatial.sh` (step 10).

## Accession table

| Cohort (label in code) | Modality | Disease(s) | Repository | Accession |
|---|---|---|---|---|
| Neyazi | single-nucleus RNA-seq + spatial | Cardiac sarcoidosis | GEO | GSE319770 / GSE319771 |
| Reichart | single-nucleus RNA-seq | Dilated & arrhythmogenic cardiomyopathy; non-failing control | EGA (via CELLxGENE) | EGAS00001006374 |
| Larson | single-nucleus RNA-seq | Hypertrophic cardiomyopathy | GEO | GSE174691 |
| Chin (2022) | single-nucleus RNA-seq | Hypertrophic cardiomyopathy (+ non-failing) | GEO | GSE181764 |
| Chin (2021) | single-nucleus RNA-seq | Non-failing control | GEO | GSE161921 |
| Liu | single-nucleus RNA-seq | Cardiac sarcoidosis; ischemic cardiomyopathy (4 CS + 3 ICM) | GEO | GSE205734 |
| Foong | spatial transcriptomics (Visium) | Cardiac sarcoidosis | GEO | GSE314910 |

## Notes

- **Open vs. controlled access.** All cohorts above are open-access and drive the reproducible
  pipeline end-to-end. The Reichart snRNA-seq matrix is retrieved from its open CELLxGENE
  distribution; the underlying EGA study (EGAS00001006374) is controlled-access and is not required
  for this pipeline.
- **Foong Visium accession.** The spatial series is GEO **GSE314910** (12 CS myocardial Visium
  samples, CS-only). The step-10 downloader takes the accession from the `FOONG_GSE` environment
  variable (defaulting to `GSE314910`) so it can be repointed without editing the script.
- **De-identification.** All samples were de-identified by the original data generators prior to
  public deposition. Consult each accession's associated publication and repository record for
  consent and data-use terms before reuse.

See the repository `README.md` for the same accession summary alongside the pipeline overview.
