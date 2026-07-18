# Cluster setup runbook: environment + smoke test

This is the **battle-tested** path: every step here reflects a real constraint commonly hit on an
HPC/SLURM cluster (login-node memory cap, shared-filesystem quota, an old system conda, the R↔Python
bridge). Follow it top to bottom.

Assumes the repo is deployed to `$PROJ` (clone or copy the bundle there) and you have edited
`env/paths.sh` to set your `LAB` (and `KID` if it differs from `$USER`). Each pipeline runner
creates the output subdirs it needs via `mkdir -p`, so no separate scaffold step is required.
Replace `<conda-module>`, `<lab>`, `<kid>` where noted.

> **TL;DR of the gotchas** (details inline below):
> 1. An old system `conda` may not solve this env and a stray `mamba` on PATH can be an unrelated test-runner: install a project **Miniforge**.
> 2. Login nodes often have a **small memory cap** (e.g. ~4 GB) → the env solve OOMs there → build inside a **32 GB `srun`**.
> 3. Persistent storage (`/gpfs/data`) often has a **space + inode quota** that a conda env can blow → put the env + package cache on **scratch** for now; raise the lab quota for the long term.
> 4. The R `anndata` reader uses **reticulate** → must point it at the env's Python (handled by `modules.sh`).

---

## 0. Paths

```bash
ssh <kid>@<cluster-login-host>
source /gpfs/data/<lab>/cs-cmsignature/env/paths.sh
echo "$PROJ $SCRATCH_ROOT $CONDA_ENV"
```

## 1. tmux (so long builds survive disconnects)

```bash
tmux new -s cs-cmsignature
source $PROJ/env/paths.sh
```
Detach: **Ctrl-b** then **d**. Reattach: `tmux attach -t cs-cmsignature`. List: `tmux ls`.

## 2. Install project-local Miniforge

If the system `conda` has no libmamba it may hang, and a `mamba` already on PATH can be an unrelated
Python test runner. Install your own:

```bash
cd $PROJ/env
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p $PROJ/env/miniforge
source $PROJ/env/miniforge/etc/profile.d/conda.sh
which mamba          # must be $PROJ/env/miniforge/bin/mamba
```

## 3. Put the env + package cache on scratch (quota), and add both to paths.sh

The persistent-storage quota (space **and** inodes) may not hold a conda env alongside large data.
Until the lab quota is raised, host the env on scratch. Edit `$PROJ/env/paths.sh` so these are the
values:

```bash
export CONDA_ENV="$SCRATCH_ROOT/conda/cs-cmsignature"     # was $PROJ/env/conda/cs-cmsignature
export CONDA_PKGS_DIRS="$SCRATCH_ROOT/conda_pkgs"   # add this line (package cache off /gpfs/data)
```

Then re-source: `source $PROJ/env/paths.sh`.

> **Scratch is purged.** This is fine for getting running now. For the real runs, ask
> your HPC support team to raise the lab `/gpfs/data` **inode** quota (a conda env is ~150k+ files),
> then rebuild the env under `$PROJ` and point `CONDA_ENV` back. Also keep any large reference data on
> scratch or its own allocation, not eating the lab data quota.

## 4. Build the env inside a 32 GB interactive job (NOT on the login node)

The solve needs more than the login node's ~4 GB cap, so it OOMs there (exit 137). Get a node with
memory and egress:

```bash
srun --partition=cpu_short --cpus-per-task=4 --mem=32G --time=03:00:00 --pty bash
source $PROJ/env/paths.sh
source $PROJ/env/miniforge/etc/profile.d/conda.sh
export CONDA_PKGS_DIRS=$SCRATCH_ROOT/conda_pkgs

# confirm this node can reach the package servers:
curl -sSI --max-time 15 https://conda.anaconda.org/bioconda/noarch/repodata.json | head -1

mamba env create -f $PROJ/env/environment.yml -p $CONDA_ENV 2>&1 | tee $PROJ/env/build.log
```

If you ever see truncated-write / "Failed to create dir" errors mid-extraction, that's a **quota**
hit: clean the partial cache (`rm -rf $CONDA_PKGS_DIRS/* $CONDA_ENV`) and make sure step 3 pointed
both at scratch, then rebuild.

## 5. Verify

```bash
conda activate $CONDA_ENV
python -c "import anndata, scanpy, scvi; print('py OK')"
Rscript -e 'library(optparse); library(anndata); library(limma); library(celda); cat("R OK\n")'
```
Both `py OK` and `R OK` must print. The R line also confirms the reticulate bridge (it imports the
Python `anndata` via the env's Python).

## 6. Confirm modules.sh (already correct in the bundle)

`env/modules.sh` sources the project Miniforge, activates `$CONDA_ENV`, and exports
`RETICULATE_PYTHON=$CONDA_PREFIX/bin/python` so the `.R` scripts reach the env's Python. If you edited a
condensed copy on the cluster earlier, re-sync the bundle's version:

```bash
rsync -av ./env/modules.sh <kid>@<cluster-login-host>:$PROJ/env/
```
Verify the Singularity line matches reality: `module avail singularity`.

## 7. Smoke test

Always submit from `$PROJ/slurm`:

```bash
cd $PROJ/slurm
sbatch smoke_test.sbatch
squeue -u $USER
cat cs_smoke_*.out | tail -20      # newest job: want PASS=11  FAIL=0
```

The smoke test points `gates`/`pseudobulk_de` at the synthetic ground-truth `cell_type` column (a
2-epoch scANVI on toy data can't produce reliable `predicted_cell_type`); the real pipeline keeps the
default `predicted_cell_type`.

Optional heavier legs:
```bash
SMOKE_DECONTX=1   sbatch smoke_test.sbatch
SMOKE_CELLBENDER=1 sbatch --partition=gpu4_short --gres=gpu:1 smoke_test.sbatch
```

---

## Troubleshooting reference (what each symptom means)

| Symptom | Cause | Fix |
|---|---|---|
| Job FAILs at `00:00:00`, ~6 MB RSS | `BASH_SOURCE` points at Slurm spool dir; `paths.sh` not sourced | scripts now anchor on `SLURM_SUBMIT_DIR`: submit from `$PROJ/slurm` |
| `mamba: unrecognized arguments: -p` | stray `mamba` test-runner on PATH | use `$PROJ/env/miniforge/bin/mamba` |
| Solve hangs at "Collecting package metadata" | system conda 4.13, classic solver | use Miniforge mamba |
| `mamba` exits 137, `oom-kill ... task=mamba` | login-node ~4 GB cap | build in a 32 GB `srun` |
| "Failed to create dir" / truncated writes mid-extract | GPFS space/inode quota | env + cache on scratch; clean partials; raise quota |
| `ModuleNotFoundError: No module named 'anndata'` (from R) | reticulate using wrong Python | `RETICULATE_PYTHON=$CONDA_PREFIX/bin/python` (in `modules.sh`) |
| `there is no package called 'optparse'` | R deps missing | `r-optparse`, `r-anndata`, `bioconductor-dropletutils` now in `environment.yml` |
| `CM nuclei: 0` / `'counts' must contain at least one value` | 2-epoch scANVI predicts no CMs on toy data | smoke uses `--celltype_col cell_type` |
| `scib metrics skipped` | scanpy/scib version drift | harmless; optional kBET/LISI only: pin scib if you want them |

## Once FAIL=0: moving to real data

1. Raise the lab `/gpfs/data` inode quota (your HPC support team), rebuild env under `$PROJ`, repoint `CONDA_ENV`.
2. Stage real data per `README.md` (GEO + Reichart CELLxGENE) and fill
   `raw/sample_meta.tsv`.
3. `build_manifest.py`, set `02_cellbender.sbatch --array=1-N%4`, then `bash submit_all.sh`.
4. Review the two gates (`$PROJ/de/gates/`) before trusting any DEG.
