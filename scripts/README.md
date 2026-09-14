# Analysis scripts

Every script in this tree is part of the manuscript pipeline. Run them with
`bash scripts/spatial_first/run_all.sh` from the repo root; see `README.md` for the
interpreter overrides the two-environment split needs.

| Script | Inputs | Output |
|---|---|---|
| `spatial_first/00_export_from_raw.py` | raw Visium objects (not redistributed) | `data/spatial_first/*` |
| `spatial_first/01_limma_q1q2.R` | pseudobulk, detection fractions, section metadata | `q1_raw_T2000_d25.tsv`, `q2_raw_T2000_d25.tsv` |
| `spatial_first/01b_limma_twogroup.R` | same, plus a two-group metadata file | any single two-group contrast (used by 05) |
| `spatial_first/02_leave_one_out.R` | same, `--mode section\|disease` | `loo/q1_raw_drop*.tsv` |
| `spatial_first/03_locus_qc.py` | Q2 results, gene spans, segdup track, annotation | `q2_locus_qc.tsv` |
| `spatial_first/04_snrna_pool.py` | Q1 results, refits, `data/snrna/snrna_CS_vs_*.tsv` | `q1_stable_down_115.tsv`, `q1_core_snrna_pergene.tsv`, `q1_core_snrna_concordance.tsv`, `snrna_gene_length_regression.tsv` |
| `spatial_first/05_controls.py` | pseudobulk + the control metadata files | `q1_control_comparison.tsv`, `q1_noMT_*.tsv`, `ctrl_panelB_vs_panelA.tsv` |
| `spatial_first/06_published_reference.py` | `data/published/*`, Q1 results | `q1_published_nf_pergene.tsv`, `q1_published_nf_summary.tsv` |
| `figure1_build.py` | `results/tables/figure1_panel{A,B}_values.tsv` | `figures/Figure1.{png,pdf}` |
| `manuscript_results.py` | everything above | `results/manuscript_results.tsv`, non-zero exit on any mismatch |

Conventions: R stages take `--pb/--det/--meta/--thresh/--out` and write
tab-separated tables with a `gene` column; Python stages default their paths to the
committed locations and need no arguments. `05_controls.py --skip-fits` rebuilds the
comparison tables without refitting, and `03_locus_qc.py --refresh` re-fetches the Ensembl
annotation instead of using the committed copy.

## Stage 00 is not run by the driver

`00_export_from_raw.py` needs the raw Visium objects, which are not redistributed (see
`docs/data_availability.md`). Its outputs are committed under `data/spatial_first/` and
`run_all.sh` starts from them. It is kept here because it documents how those inputs were
derived: spot-level cardiomyocyte/immune classification, the 0.5 mm exclusion radius around
immune-enriched spots, the section-evaluability rule, and the per-section pseudobulk sum.

