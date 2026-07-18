# Cardiac Sarcoidosis Cardiomyocyte Signature

Integrative single-nucleus and spatial transcriptomic datasets were analyzed to suggest a
**cardiomyocyte transcriptional signature associated with cardiac sarcoidosis (CS)**. Candidate
genes were nominated by comparing CS to multiple cardiomyopathies and non-failing hearts with
publicly available snRNA-seq. The genes were validated in two independent datasets: a
CS-versus-ICM single-nucleus comparison (Liu et al., Circulation Research, GSE205734) and spatial
(Visium) transcriptomics of intact CS myocardium (Foong et al., Journal of Cardiac Failure:
Intersections, GSE314910). A within-patient mixed-effects distance-to-lesion analysis established
the cell-autonomous nature of the gene signature. The pipeline runs on an HPC cluster with SLURM
job management system via conda and singularity.

This repository of analysis code accompanies the Research Letter *"Toward a Molecular Diagnosis of
Cardiac Sarcoidosis with a Cardiomyocyte-Intrinsic Transcriptional Signature"*. Analysis code and derived result tables are included.

## Key result

The analysis suggests a four-gene cell-autonomous cardiomyocyte signature of cardiac sarcoidosis:

> **GJB7, TNNI3K, MLIP, PANK1**

The validation design is tiered and aware of confounds:

- The **cross-cardiomyopathy discovery contrasts** (CS versus DCM / ARVC / HCM / non-failing,
  assembled across cohorts) are **cohort-confounded** as cardiac sarcoidosis is contributed by a
  single cohort in the discovery phase. These contrasts are **hypothesis-generating** and do not
  on their own support the signature claim.
- Instead, the claim is validated by two **independent tiers**:
  1. **Within-study Liu CS-vs-ICM** (GSE205734): CS and ICM are compared inside a single cohort,
     removing the cross-cohort batch confound.
  2. **Foong spatial (Visium)** (GSE314910): the signature is further validated in intact CS
     myocardium, including a within-patient distance-to-lesion mixed-effects model that establishes
     cell-autonomy (the signal tracks cardiomyocytes rather than infiltrating immune/granuloma
     content).

## Repository structure

```
/
├── README.md                # this file
├── LICENSE                  # MIT
├── .gitignore
├── .gitattributes           # Git LFS patterns for large result tables
├── environment.yml          # conda environment (also in env/)
├── env/
│   ├── environment.yml
│   ├── cellbender.def       # Singularity definition for CellBender
│   ├── paths.sh             # central path/identity config sourced by every sbatch (edit LAB)
│   └── modules.sh           # Lmod modules + conda activation sourced by every sbatch
├── scripts/                 # all analysis scripts (Python + R) + scripts/README.md
├── slurm/                   # SLURM runners, steps 00–11 + smoke_test + submit_all.sh
├── docs/
│   ├── pipeline.md          # pipeline overview
│   ├── setup_cluster.md     # cluster setup notes
│   ├── next_steps.md
│   └── data_availability.md # accession table
├── metadata/
│   ├── neyazi_cs_zones_by_patient.tsv
│   ├── neyazi_cs_zones_by_section.tsv
│   └── sample_meta.TEMPLATE.tsv
├── figures/
│   ├── Figure1_CS_cardiomyocyte.png
│   └── Figure1_CS_cardiomyocyte.pdf
└── results/                 # derived result DATA TABLES (large files via Git LFS)
```

## Data availability

All datasets are publicly available and de-identified. No raw data is redistributed in this
repository; only analysis code and derived result tables are included.

| Cohort (label in code) | Modality | Disease(s) | Accession | Source publication |
|---|---|---|---|---|
| Neyazi | single-nucleus RNA-seq + spatial | Cardiac sarcoidosis | GSE319770 / GSE319771 | Neyazi et al. (in press) |
| Reichart | single-nucleus RNA-seq | Dilated & arrhythmogenic cardiomyopathy; non-failing control | EGAS00001006374 | Reichart et al., *Science* 2022, [10.1126/science.abo1984](https://doi.org/10.1126/science.abo1984) |
| Larson | single-nucleus RNA-seq | Hypertrophic cardiomyopathy | GSE174691 | Larson et al., *Sci Rep* 2022, [10.1038/s41598-022-08561-x](https://doi.org/10.1038/s41598-022-08561-x) |
| Chin (2022) | single-nucleus RNA-seq | Hypertrophic cardiomyopathy (+ non-failing) | GSE181764 | Codden et al., *Int J Mol Sci* 2022, [10.3390/ijms23020946](https://doi.org/10.3390/ijms23020946) |
| Chin (2021) | single-nucleus RNA-seq | Non-failing control | GSE161921 | Larson et al., *BMC Med Genomics* 2021, [10.1186/s12920-021-01011-z](https://doi.org/10.1186/s12920-021-01011-z) |
| Liu | single-nucleus RNA-seq | Cardiac sarcoidosis; ischemic cardiomyopathy (4 CS + 3 ICM) | GSE205734 | Liu et al., *Circ Res* 2022, [10.1161/CIRCRESAHA.121.320449](https://doi.org/10.1161/CIRCRESAHA.121.320449) |
| Foong | spatial transcriptomics (Visium) | Cardiac sarcoidosis | GSE314910 | Foong et al., *J Card Fail Intersect* 2026, [10.1016/j.yjcafi.2025.12.009](https://doi.org/10.1016/j.yjcafi.2025.12.009) |

## Pipeline / reproducibility

The pipeline is a SLURM dependency chain, run in order 00 → 11. Each step is an
`sbatch` runner in `slurm/` that calls one or more scripts in `scripts/`. The table below is
derived directly from the `slurm/*.sbatch` files.

| Step | Runner | Script(s) invoked | Purpose |
|---|---|---|---|
| 00 | `00_setup_env.sbatch` | (conda env create; `singularity pull`) | Provision the conda env and pull the CellBender Singularity image |
| 01 | `01_download.sbatch` | `download_geo.py`, `download_reichart_cxg.py`, `build_manifest.py` | Download public GEO / CELLxGENE datasets and build the sample manifest |
| 02 | `02_cellbender.sbatch` | CellBender `remove-background` (via Singularity; job array) | Ambient-RNA removal on raw count matrices |
| 03 | `03_decontx.sbatch` | `run_decontx.R` | DecontX ambient-contamination correction |
| 04 | `04_qc.sbatch` | `run_qc.py` | Per-dataset QC filtering; merge sample metadata onto obs |
| 05 | `05_scanvi.sbatch` | `run_scanvi.py` | scANVI integration across cohorts (batch correction) |
| 05b | `05b_relabel.sbatch` | `relabel_genes.py` | Relabel integrated + liu_qc var_names from Ensembl IDs to HGNC symbols (required by all symbol-space downstream stages) |
| 06 | `06_validate_gates.sbatch` | `run_gates.py` | Cell-type gate validation and QC UMAPs |
| 07 | `07_pseudobulk_de.sbatch` | `subtype_tnk.py`, `pseudobulk_cm.py`, `run_pseudobulk_de.R` | Pseudobulk aggregation and cross-cardiomyopathy differential expression (limma-voom / pyDESeq2) |
| 08 | `08_liu_standalone.sbatch` | `run_liu_standalone.py` | Within-study Liu CS-vs-ICM differential expression |
| 09 | `09_gsea_figures.sbatch` | `run_gsea_figures.py`, `annotate_foong_validation.py`, `panel_robustness.py` | GSEA, volcano/enrichment figures, Foong validation annotation, panel robustness |
| 10 | `10_foong_spatial.sbatch` | `download_foong_spatial.sh`, `foong_spatial_figures.py` | Download Foong Visium data and generate spatial figures |
| 11 | `11_cm_signature.sbatch` | `cm_disease_distance.py`, `cm_spatial_crossdisease.py`, `cm_content_normalized.py`, `foong_regional.py` | Cardiomyocyte signature: within-patient distance-to-lesion, cross-disease spatial, CM-content normalization, regional analysis |

**Smoke test.** `slurm/smoke_test.sbatch` runs the whole chain (make_synthetic_data →
build_manifest → decontx → cellbender → qc → scanvi → relabel → gates → pseudobulk_de → liu_standalone →
gsea_figures) on small synthetic data to verify wiring before committing
cluster resources.

**Full submission.** `slurm/submit_all.sh` launches steps 01–11 as an `afterok` dependency chain
(07/08 fan out in parallel after the gates step, 09 joins them, then 10 spatial and 11 CM-signature
run in sequence). Edit the `02_cellbender.sbatch` `--array` range to match your manifest sample count
before launching.

Helper scripts not wired into a numbered step are run manually: `composite_figure.py` (assembles
the final Figure 1), `combine_contrasts.py`, `diagnose_neyazi_liu.py`, `foong_panelF_candidates.py`,
`relabel_genes.py`, `make_sample_meta.py`. `_common.py` is a shared import
(logging, gene panels). See `scripts/README.md` for the full per-script contract.

## Environment setup

The analysis environment is a conda environment defined in `environment.yml`:

```bash
conda env create -f environment.yml
conda activate cs-cmsignature
```

CellBender runs from a Singularity image built from `env/cellbender.def` (or pulled directly, as in
step 00). The pipeline targets a SLURM cluster.

`env/paths.sh` centralizes all paths and cluster identity; **edit the `LAB` variable** (your
`/gpfs/data/<lab>` group) before running. `env/modules.sh` handles Lmod module loading and conda
activation. Both are sourced by every `sbatch` runner. The `#SBATCH --account` / `--partition`
lines and any cluster-specific paths are placeholders (`YOUR_ACCOUNT`, `YOUR_PARTITION`): set them
for your site (verify partitions with `sinfo -s`).

## How to reproduce from public data

1. **Download** (step 01): fetch the public datasets in the accession table and build the manifest.
2. **Decontaminate** (steps 02–03): CellBender + DecontX ambient-RNA correction.
3. **QC & integrate** (steps 04–05): per-dataset QC, then scANVI integration.
4. **Gates** (step 06): cell-type gate validation.
5. **Pseudobulk DE** (step 07): cross-cardiomyopathy contrasts (hypothesis-generating).
6. **Liu within-study** (step 08): batch-clean CS-vs-ICM validation tier.
7. **GSEA / validation / robustness** (step 09): enrichment, Foong annotation, panel robustness.
8. **Spatial** (step 10): Foong Visium figures.
9. **CM signature** (step 11): within-patient distance-to-lesion, cross-disease spatial, CM-content
   normalization, regional analysis. Assemble Figure 1 with `composite_figure.py`.

## Known limitations

These are inherent to the design and data, and are documented here for anyone reusing the code:

- The **cross-cardiomyopathy discovery contrasts are cohort-confounded** (cardiac sarcoidosis from a
  single cohort → disease aliased with batch → rank-deficient pseudobulk design). They are
  hypothesis-generating; the independent Liu within-study and Foong spatial tiers carry the claim.
- **Ambient correction is non-uniform across cohorts** (different chemistries and available raw
  matrices lead to CellBender vs DecontX handling per dataset).
- `make_sample_meta.py` infers disease by **substring matching**, which has a `nicm` ⊃ `icm`
  edge case. It is safe for the datasets used here, but should be checked before reuse on other
  cohorts.

## Citation

If you use this software or its results, please cite the associated article. The manuscript DOI
will be added on publication.

## License

MIT: see [`LICENSE`](LICENSE).
