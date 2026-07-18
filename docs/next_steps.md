# NEXT_STEPS: from green smoke test to real results

Checklist for moving the validated pipeline onto real data. Tick top to bottom. Commands assume
`source $PROJ/env/paths.sh` has been run.

---

## Phase 1: Make the environment durable (do first)
The env currently lives on **scratch** (purged periodically). Move it to persistent storage before a
multi-day real run.

- [ ] Email your HPC support team: raise lab `/gpfs/data/<lab>` **inode** quota (~150k+ files
      for a conda env) and confirm ~1 TB space headroom (working set + 80 GB atlas).
- [ ] After quota raised: rebuild env under `$PROJ` (32 GB `srun` + Miniforge mamba; see SETUP §4):
      `mamba env create -f $PROJ/env/environment.yml -p $PROJ/env/conda/cs-cmsignature`
- [ ] Set `CONDA_ENV="$PROJ/env/conda/cs-cmsignature"` back in `env/paths.sh`; re-source.
- [ ] Re-verify: `python -c "import anndata,scanpy,scvi"` and `Rscript -e 'library(celda);library(limma)'`.
- [ ] (If quota is slow) you may run the real pipeline from the scratch env meanwhile: just don't let
      it idle long enough to be purged.

## Phase 2: Phase 0 access decisions (governance gate)
- [ ] Open `env/access_status.md`; keep Simonson-ICM / ARVC / Chaffin-HCM as **`atlas-transfer`**
      (default: uses the 80 GB atlas you already have; no controlled raw needed).
- [ ] Only if a reviewer requires raw co-processing: flip an arm to `raw-reprocess` and start the dbGaP
      **phs001539** DAR + DUA + restricted-dir request now (approvals take weeks).

## Phase 3: Stage real data + metadata (biggest manual lift)
- [ ] Write `raw/accessions.tsv` (dataset⇥gse⇥supp_url?) for the GEO pulls.
- [ ] Download (data-mover node or `01_download.sbatch`):
  - [ ] Neyazi CS: GSE319770 (41 raw_feature_bc_matrix H5)
  - [ ] Larson HCM: GSE174691
  - [ ] Chin2022: GSE181764 ; Chin2021: GSE161921
  - [ ] Reichart "All cells": `https://datasets.cellxgene.cziscience.com/6d0b8b8a-2b22-4578-b8ee-8ca4ab893594.h5ad`
  - [ ] Liu: GSE205734 (5′ arm)
- [ ] Lay out exactly as the scripts expect:
      `raw/<dataset>/<sample>/raw_feature_bc_matrix/{matrix.mtx.gz,features.tsv.gz,barcodes.tsv.gz}`
      for the GEO sets; `raw/reichart/reichart_cxg.h5ad`; `raw/liu/...`.
- [ ] **Fill `raw/sample_meta.tsv`** (template provided: `raw/sample_meta.TEMPLATE.tsv`). One row per
      sample with `disease, procurement, region`. This is the piece only you have, and the DE model +
      procurement gate depend on it. Get Neyazi's 41 CS samples and the procurement labels right.
- [ ] Verify Reichart provenance note (raw counts in `raw.X`, not CellBender: already confirmed).
- [x] **`sample_meta.tsv` → obs merge (DONE).** `run_qc.py` now merges disease/anatomy/sex/individual
      onto each cell by `sample_id` (which `run_decontx.R` writes into obs), non-destructively (Reichart's
      own obs is left intact). `04_qc.sbatch` passes `--sample_meta $RAW/sample_meta.tsv`. Just make sure
      the file is filled before launch.
- [x] **Disease-label harmonization (DONE).** `run_qc.py` maps every dataset's `disease` to one vocab
      (`CS/DCM/ICM/ARVC/HCM/NF/NCC`) via `_common.normalize_disease`: Reichart's CELLxGENE ontology
      strings ("dilated cardiomyopathy" → `DCM`, "normal" → `NF`, etc.) and the GEO short codes all land
      in the same space. Unrecognized labels pass through with a warning; add them to
      `_common.DISEASE_ALIASES` if any show up in the logs.

## Phase 4: Manifest, array sizing, launch
- [ ] `python scripts/build_manifest.py --raw $PROJ/raw --out $PROJ/raw/manifest.tsv`
- [ ] Sanity-check the knee `expected_cells` look reasonable on the real barcode-rank curves.
- [ ] Set `slurm/02_cellbender.sbatch --array=1-N%4` where N = `$(( $(wc -l < $PROJ/raw/manifest.tsv) - 1 ))`.
- [ ] Reconcile partitions/walltimes once more on the live cluster (`sinfo -s`).
- [ ] Launch: `cd $PROJ/slurm && bash submit_all.sh`
- [ ] Monitor: `squeue -u $USER` ; `sacct -j <jid> --format=JobID,State,Elapsed,MaxRSS`.
- [ ] Watch **CellBender** (GPU array, ~13–26 h: the wall-clock driver) and **scANVI** (OOM risk at
      64 GB on ~770k nuclei → resubmit to an `fn_*` high-mem node, `--mem=192G`).

## Phase 5: Review the two gates (the scientific payoff)
- [ ] **Decontamination QC gate** (`de/gates/decontam_qc_gate.tsv`): lineage-ambient panel ≈ 0 in every
      contrast; CS-vs-DCM drops from the old +1.86 toward 0.
- [ ] **Dual-count concordance**: wire raw CellRanger counts into `run_gates.py` (currently a TODO
      stub) and require DEGs concordant in both count versions + background heuristic < 0.5.
- [ ] **Procurement-robustness gate**: candidates must survive across DCM, Simonson-ICM, ARVC: not just
      the mismatched HCM-myectomy / NF. Expect the short real signature (TNNI3K, GJB7); CM-dedifferentiation
      panel stays a labeled non-specific control.
- [ ] **Liu co-primary**: check overlap (hypergeometric) + direction concordance with the 3′ primaries.
- [ ] **Cross-disease projection**: where CS lands vs DCM/HCM/ICM/ARVC; COL22A1⁺ (DCM) vs TNC⁺ (ICM)
      fibroblast niche occupancy.
- [ ] Only report DEGs that pass **both** gates.

## Loose ends (whenever)
- [ ] Dual-count concordance raw-counts wiring (`run_gates.py`).
- [ ] Decide whether to pin `scib` for kBET/LISI integration metrics (optional; gates work without it).
- [ ] After quota raise, move the 80 GB atlas off any quota-limited space if needed.
