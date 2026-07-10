# Plan v3.0 (BigPurple) — Unified scANVI Pipeline: Cardiac Sarcoidosis vs. Mimics

**Target platform: NYU Langone BigPurple HPC.** This revises Plan v2.1 for execution entirely on
BigPurple. **The science is unchanged** — same datasets, same uniform-decontamination logic, same two
gates (decontamination QC + procurement robustness), same scANVI integration and HeartMap-model DE.
What changes is the *infrastructure*: Slurm scheduling, Lmod modules, Singularity containers (no
Docker on BigPurple), GPFS storage layout, GPU-partition CellBender, and controlled-access data
staging under NYU Langone's data-governance process.

If you want the underlying scientific rationale (the procurement confound, the TNNI3K/GJB7 correction,
the HeartMap incorporation), see Plan v2.1 — it is carried forward verbatim in intent. This document
is the **execution layer**.

---

## 0. What BigPurple changes (orientation)

| Concern | bionmi (old) | BigPurple (new) |
|---|---|---|
| Scheduler | ad hoc / local | **Slurm** (`sbatch`, `srun`, `salloc`, job arrays, dependencies) |
| Software | local installs | **Lmod modules** (`module load`) + **conda/mamba** env + **Singularity** for CellBender |
| Containers | Docker | **Singularity/Apptainer only** (build/pull on a build node or pull `docker://`) |
| GPU | single/none | **gpu4 / gpu8 partitions** (V100 32 GB / A100 where available), Slurm `--gres=gpu:N` |
| Storage | local disk | **GPFS**: `/gpfs/data/<lab>` (persistent), `/gpfs/scratch/<kid>` (fast, purged), `/gpfs/home/<kid>` |
| Parallelism | sequential | **job arrays** (`--array`) — one CellBender task per sample, 4–8 concurrent |
| Controlled data | — | **dbGaP phs001539 / EGA** staged only after DUA + NYU data-governance sign-off |

**Three rules that drive every choice below:**

1. **Heavy I/O and intermediates live on `/gpfs/scratch`** (flash, fast) but scratch is **purged**
   (typically inactive-file purge ~weeks) — so **persist every durable artifact** (CellBender h5,
   DecontX counts, the integrated `.h5ad`, DE tables) **back to `/gpfs/data/<lab>`** as each stage
   finishes. The original v2.1 "persist to `/mnt/shared-workspace/`" maps to **`/gpfs/data/<lab>/heartmap-cs/`**.
2. **CellBender runs in a Singularity container on a GPU partition**, one sample per array task.
3. **Controlled raw data (Chaffin/Simonson dbGaP, Reichart EGA) cannot be downloaded until access is
   granted** — Phase 0 below makes that an explicit gate, and the pipeline is designed to run the
   **open-data arms first** and fold controlled arms in later without rework.

> **Verify cluster specifics before submitting.** Partition names and walltime caps change. Run
> `sinfo -s`, `sacctmgr show qos`, and `module avail` on a login node and reconcile against the
> placeholders in this plan. Where this doc says e.g. `gpu4_medium`, confirm it exists and check its
> `MaxTime`. HPC support: **your HPC support team** / the HPC Portal.

---

## 1. Storage layout (GPFS)

Set these once in `env/paths.sh` and source everywhere. Replace `<lab>` and `<kid>` (your KerberosID).

```
PROJ=/gpfs/data/<lab>/heartmap-cs           # persistent project root (deliverables live here)
SCRATCH=/gpfs/scratch/<kid>/heartmap-cs     # fast, PURGED — intermediates only
HOME_ENV=/gpfs/home/<kid>                    # conda/mamba + Singularity cache

$PROJ/
  raw/            # downloaded matrices (open data); controlled data → $PROJ/controlled (see Phase 0)
  controlled/     # dbGaP phs001539 + EGA — restricted ACLs, governance-approved location
  cellbender/     # CellBender h5 outputs (persisted from scratch)
  decontx/        # DecontX-corrected integer counts (persisted)
  qc/             # per-dataset QC'd h5ad
  integrated/     # scANVI model + integrated h5ad
  de/             # pseudobulk + DE tables (per contrast)
  figures/        # final figures
  logs/           # Slurm .out/.err (or keep in $SCRATCH/logs and copy)
  env/            # environment.yml, modules.sh, paths.sh, *.sif build recipes
  slurm/          # all sbatch scripts

$SCRATCH/
  cellbender_work/  decontx_work/  scanvi_work/   # transient compute scratch
```

**Quota check before you start:** `gpfsquota` / `df -h /gpfs/data/<lab>` and confirm headroom.
Budget ~**0.5–1 TB** working set (52+ raw matrices, CellBender h5s, integrated ~770k-nucleus h5ad,
plus the HeartMap atlas). Keep the 3.7 TB Reichart EGA BAM *out of scope* — we use the CELLxGENE
matrix, exactly as in v2.1.

---

## 2. Software environment

Three layers, in priority order. Build/setup once via `slurm/00_setup_env.sbatch`.

**(a) Conda/mamba env (`env/environment.yml`)** — the Python/R analysis stack: scanpy, scvi-tools,
scArches, anndata, cellxgene-census, pydeseq2, gseapy, scDblFinder, scib, scCODA; R limma/edgeR,
celda(DecontX), DESeq2. BigPurple provides `miniconda3`/`anaconda3` via Lmod; create the env on
`/gpfs/home` or `/gpfs/data` (not scratch — it'd get purged).

```bash
module load miniconda3            # confirm exact name: module avail conda
mamba env create -f env/environment.yml -p $PROJ/env/conda/heartmap
```

**(b) Singularity image for CellBender (`env/cellbender.def` or a pulled `.sif`)** — CellBender +
CUDA. Pull on a node with internet (login/data-mover) into the Singularity cache, then run with
`--nv` for GPU passthrough:

```bash
module load singularity/3.x
singularity pull $PROJ/env/cellbender.sif docker://us.gcr.io/broad-dsde-methods/cellbender:0.3.0
# GPU run pattern inside the array job:
singularity exec --nv -B $SCRATCH,$PROJ $PROJ/env/cellbender.sif cellbender remove-background ...
```

> Why container for CellBender and conda for the rest: CellBender's CUDA/torch pinning fights the
> scvi-tools stack. Isolating it in a `.sif` avoids dependency conflicts and gives a reproducible GPU
> runtime. scvi-tools scANVI training runs fine from the conda env (GPU optional — see §6).

**(c) Lmod modules (`env/modules.sh`)** — `singularity`, `miniconda3`, and (if you build R outside
conda) `r/4.x`. Keep module loads in one sourced file so every sbatch script is identical.

---

## 3. Phase 0 — Access audit & data governance (GATE, do this first)

This is the load-bearing change from v2.1. On BigPurple, controlled data has a *process*, not just a
download.

**3.1 Open data (no gate) — start here, pipeline runs end-to-end on these:**
- Neyazi CS (GSE319770), Larson HCM (GSE174691), Chin2022 (GSE181764), Chin2021 (GSE161921) — GEO, open.
- Reichart (CELLxGENE matrix) — open. **Confirmed:** collection `e75342a8-0f3b-4ec5-8ee1-245a23e0f7cb`
  (Reichart 2022, *Science* abo1984; EGA EGAS00001006374). Use the **"All cells"** dataset (881,081
  nuclei, ~6.4 GB): `https://datasets.cellxgene.cziscience.com/6d0b8b8a-2b22-4578-b8ee-8ca4ab893594.h5ad`.
  Provenance answered — raw counts in `raw.X`, `.X` normalized, **not** CellBender; run uniform DecontX
  on `raw.X`. Chemistry is **mixed 3′ v2 + v3** (audit table assumed pure v3); lymphocytes lumped.
- HeartMap **integrated atlas** (Broad Single Cell Portal **SCP3689**) — check if openly downloadable
  (atlases usually are even when per-study raw is gated). If yes, this gives you the **label/signature
  transfer** path for ICM/ARVC/LV-HCM **without** touching controlled raw.

**3.2 Controlled raw (gated) — Chaffin & Simonson via dbGaP phs001539; Reichart raw via EGA:**
You indicated you need access guidance, so the path on BigPurple is:

1. **Confirm whether you even need controlled raw.** If SCP3689 (the harmonized atlas) is
   downloadable, v2.1's *default* is **transfer/project, not raw co-embed** — so you may not need
   phs001539 at all. Decide per arm. Recommended default: **atlas-only**, fold in raw later only if a
   reviewer demands it.
2. **dbGaP phs001539 (controlled):** requires an approved **Data Access Request (DAR)** through your
   PI's eRA Commons, with an institutional signing official. Approval → you receive a **dbGaP
   repository key (`.ngc`)**.
3. **NYU Langone data governance:** controlled-access human genomic data on BigPurple must land in a
   **governance-approved, access-restricted location** with locked-down GPFS ACLs and an approved Data
   Use Agreement on file. **Email your HPC support team before staging** to set up the restricted
   directory (`$PROJ/controlled`) and confirm dbGaP/SRA-toolkit use is permitted on the data-mover
   nodes. Do not place controlled data in world/group-readable space.
4. **Download mechanics (only after 1–3):** SRA Toolkit `prefetch`/`fasterq-dump` with the `.ngc` key
   for dbGaP; **`pyega3`** for EGA. Run on **data-mover nodes**, not login nodes. (Reichart EGA raw
   stays out of scope per v2.1 — BAM-only, 3.7 TB; we use the CELLxGENE matrix.)

**Phase 0 decision recorded in `env/access_status.md`:** for each of {Simonson ICM, ARVC, Chaffin
LV-HCM} mark `atlas-transfer` (default) or `raw-reprocess` (requires phs001539). The pipeline branches
on this file.

---

## 4. Phase 1 — Dataset selection & the procurement confound

**Unchanged from v2.1.** Chemistry/procurement audit, the procurement × disease confound table, and
inclusion (integrated 3′ scANVI arm + standalone Liu 5′ co-primary + HeartMap reference option) all
carry over exactly. The only BigPurple-relevant note: build a **single sample manifest TSV**
(`$PROJ/raw/manifest.tsv`) listing every sample with columns `sample_id, dataset, accession, chemistry,
procurement, region, raw_path, expected_cells_seed`. **Every Slurm array job indexes into this
manifest** — this is what makes the array sizes self-documenting and removes the "52 vs 41" ambiguity
(the array length = number of rows with raw droplets, computed at submit time, not hard-coded).

---

## 5. Phase 2 — Unified processing pipeline (Slurm form)

Each step below maps to one script in `slurm/`. The scientific content (uniform decontamination, QC
thresholds, harmonized annotation, feature intersection, the two gates) is **identical to v2.1**; only
the execution wrapper is new.

### Step 1 — Acquisition + Reichart verification (`01_download.sbatch`)
- Partition: `data_mover` / `cpu_short`. CPU, low mem, network egress.
- GEO pulls (Neyazi, Larson, Chin×2) + Reichart CELLxGENE matrix → `$PROJ/raw/`.
- **Verify Reichart decontamination state** (inspect provenance; do not assume CellBender). Record in
  `$PROJ/raw/reichart_provenance.md`.
- Controlled arms: only if `access_status.md` says `raw-reprocess` (data-mover, `$PROJ/controlled`).

### Step 2a — CellBender from raw (`02_cellbender.sbatch`, **GPU job array**)
The four GEO datasets (raw droplets available). **One array task per sample**, identical params.

```bash
#SBATCH --partition=gpu4_medium          # verify name + MaxTime with sinfo
#SBATCH --gres=gpu:1
#SBATCH --array=1-N%4                     # N = #samples w/ raw droplets; %4 = 4 concurrent
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
```
- `N` and per-sample `--expected-cells` / `--total-droplets-included` are read from `manifest.tsv`
  (knee estimate computed programmatically — **no hand-tuning**, same algorithm everywhere).
- Runs `singularity exec --nv $PROJ/env/cellbender.sif cellbender remove-background ...` with the
  exact v2.1 params (`--fpr 0.01 --epochs 150 --total-droplets-included 25000` seed, `--cuda`).
- Write to `$SCRATCH/cellbender_work/<sample>/`, then **rsync the `.h5` to `$PROJ/cellbender/`** on success.
- **Concurrency `%4`** respects GPU availability; raise to `%8` if you have gpu8 headroom.

### Step 2b — Uniform DecontX pass (`03_decontx.sbatch`, CPU array, ALL FIVE datasets)
The consistency fix: identical DecontX (celda) on the **filtered** matrix of every dataset (CellBender
output for the 4 GEO; CELLxGENE filtered matrix for Reichart). Round to integers for pseudobulk.
- Partition `cpu_medium`, `--mem=16G`, `--array` over datasets. Outputs → `$PROJ/decontx/`.

### Step 2c — QC + doublets (`04_qc.sbatch`)
v2.1 thresholds verbatim: min_genes ≥200, max_genes ≤8000, MT% <5%, CellBender cell-prob >0.5;
scDblFinder/scrublet on decontaminated counts; genes in ≥3 nuclei/sample; protein_coding + lncRNA.
`cpu_medium`, `--mem=32G`, per-dataset array → `$PROJ/qc/`.

### Step 3 — Harmonized annotation
Reichart CELLxGENE Cell-Ontology labels as reference; scANVI semi-supervised transfer. **B-cell split
from T/NK** (v2.1 change retained). Folded into the scANVI step.

### Step 4 — Gene-ID harmonization & feature intersection
Map to HGNC (GRCh38), resolve dup suffixes, intersect across the five 3′ datasets (~15–18k genes).
Light CPU job or run inline at the head of Step 5.

### Step 5 — scANVI integration (`05_scanvi.sbatch`)
`.X` = post-DecontX counts; `batch_key=study`, `labels_key=cell_type_harmonized`,
`unlabeled_category="Unknown"`. Liu (5′) **not** co-embedded.
```bash
#SBATCH --partition=gpu4_short    # scANVI trains much faster on GPU; CPU fallback = cpu_long
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G                 # ~770k nuclei; bump to a fn_* high-mem node if OOM
#SBATCH --time=04:00:00
```
- `model.train(max_epochs=400, early_stopping=True)`; save model + latent + predictions to
  `$PROJ/integrated/`. (scvi-tools uses GPU automatically if visible; v2.1 said "CPU ok" — on
  BigPurple, GPU is faster and you have access, so default GPU with CPU fallback.)

### Step 6 — Validation + TWO gates (`06_validate_gates.sbatch`)
- **6a** integration validation (UMAP, kBET/LISI). **6b** decontamination QC gate incl. **dual-count
  concordance** (CellBender/DecontX vs. CellRanger raw) + background heuristic <0.5. **6c**
  procurement-robustness gate across explant-matched mimics. **Corrected gene expectations
  (TNNI3K/GJB7; CM-dedifferentiation = technical control only)** carry over verbatim. `cpu_medium`,
  `--mem=32G`.

### Step 7 — CM subset & pseudobulk DE (`07_pseudobulk_de.sbatch`)
Extract CMs; pseudobulk per patient (sum → integer). **HeartMap DE model**
`Expr ~ 0 + disease + study + anatomy + sex` with limma-voom + `duplicateCorrelation(individual)`.
Primary contrasts: CS-vs-DCM, CS-vs-ICM(Simonson 3′), CS-vs-ARVC; CS-vs-ICM(Liu) co-primary (Step 7b).
apeglm shrinkage + dupCor + ≥75% concordance + procurement gate. IEG/HSP + HeartMap-DCM generic-HF
panels. `cpu_medium`, `--mem=16G` (R + Python).

### Step 7b — Liu CS-vs-ICM standalone co-primary (`08_liu_standalone.sbatch`)
5′ arm, **not co-embedded**, recomputed fresh through QC+DecontX+pseudobulk+apeglm. Small CPU job.

### Step 8 — GSEA + figures (`10_gsea_figures.sbatch`)
Pre-ranked gseapy (Hallmark, GO BP, KEGG, Reactome) on procurement-matched primaries + high-confidence
set. All v2.1 figures/tables incl. decontamination + procurement-robustness + HeartMap-positioning.
`cpu_short`.

---

## 6. Phase 3 — Compute strategy mapped to BigPurple

| Step | Script | Partition (verify) | GPU | Mem | Walltime | Parallelism |
|---|---|---|---|---|---|---|
| Download + Reichart verify | 01 | data_mover / cpu_short | — | 8 G | ~1 h | serial |
| CellBender (per-sample) | 02 | **gpu4_medium** | 1 | 32 G | ~2–8 h/task | `--array=1-N%4` |
| DecontX (all 5) | 03 | cpu_medium | — | 16 G | ~15 min/ds | array over datasets |
| QC + doublets | 04 | cpu_medium | — | 32 G | ~30 min/ds | array over datasets |
| scANVI integrate | 05 | **gpu4_short** (CPU: cpu_long) | 1 | 64 G | 1–2 h | single |
| Validation + gates | 06 | cpu_medium | — | 32 G | ~45 min | single |
| Pseudobulk DE | 07 | cpu_medium | — | 16 G | ~30 min | single |
| Liu standalone | 08 | cpu_short | — | 8 G | ~15 min | single |
| GSEA + figures | 10 | cpu_short | — | 8 G | ~30 min | single |

**Wall-clock driver is still CellBender.** With `--array=...%4` (4 concurrent GPU tasks) over ~52
samples at ~2 h each ≈ **26 h of GPU array time**; with `%8` on gpu8 ≈ **13 h**. Everything downstream
is hours, not days. Total **~1–1.5 days** once data is staged and access is cleared.

**High-mem fallback:** if scANVI or concatenation OOMs at 64 G, resubmit to a **`fn_*` high-memory
partition** (BigPurple has 8 high-mem nodes) rather than chasing memory on a standard node.

---

## 7. Phase 4 — Orchestration (Slurm dependency chain)

`slurm/submit_all.sh` chains jobs with `--dependency=afterok` so the whole pipeline launches with one
command and each stage waits for clean completion of the prior:

```bash
jid1=$(sbatch --parsable 01_download.sbatch)
jid2=$(sbatch --parsable --dependency=afterok:$jid1 02_cellbender.sbatch)   # array
jid3=$(sbatch --parsable --dependency=afterok:$jid2 03_decontx.sbatch)
jid4=$(sbatch --parsable --dependency=afterok:$jid3 04_qc.sbatch)
jid5=$(sbatch --parsable --dependency=afterok:$jid4 05_scanvi.sbatch)
jid6=$(sbatch --parsable --dependency=afterok:$jid5 06_validate_gates.sbatch)
jid7=$(sbatch --parsable --dependency=afterok:$jid6 07_pseudobulk_de.sbatch)
jid8=$(sbatch --parsable --dependency=afterok:$jid6 08_liu_standalone.sbatch)   # parallel to 07
sbatch --dependency=afterok:$jid7:$jid8 10_gsea_figures.sbatch
```

**Run open-data arms first.** If controlled raw isn't cleared, the chain runs fully on the open arms;
fold Simonson/Chaffin raw in later by re-running from Step 2 for those samples only (the manifest +
array design makes this incremental). Default per `access_status.md` is **atlas-transfer**, which needs
no controlled raw at all.

---

## 8. Caveats — BigPurple-specific (in addition to all v2.1 caveats)

1. **All v2.1 scientific caveats stand** (Reichart can't be CellBender-re-run from raw; procurement is
   structurally confounded; HeartMap has no CS; single-CS-cohort limit; HeartMap lumps lymphocytes;
   Neyazi GEO has no published labels). None of these are fixed by changing clusters.
2. **Scratch purge risk.** The single biggest BigPurple-specific failure mode: losing CellBender/DecontX
   intermediates to a scratch purge mid-project. **Every script persists durable outputs to `$PROJ`
   immediately on success** — do not leave deliverables only on `/gpfs/scratch`.
3. **Controlled-data governance is a hard gate, not a formality.** dbGaP phs001539 on BigPurple needs an
   approved DAR + DUA + a governance-approved restricted directory. Coordinate with
   your HPC support team *before* downloading. The plan is built so you never block on this — open
   arms run first; atlas-transfer is the default.
4. **Partition names/limits drift.** Treat every `gpu4_*`/`cpu_*`/`fn_*` name and walltime here as a
   placeholder to reconcile against `sinfo -s` / `sacctmgr show qos` on the live cluster.
5. **Singularity, not Docker.** CellBender's Broad image is pulled as a `.sif`; `--nv` is required for
   GPU passthrough; bind-mount `$SCRATCH` and `$PROJ` with `-B`.

---

## 9. Execution order (BigPurple)

1. **Phase 0 access audit** — fill `env/access_status.md` (atlas-transfer default); email
   your HPC support team if any arm is `raw-reprocess`; verify SCP3689 downloadability.
2. **Provision env** — `module load` set; `mamba env create` on `$PROJ/env`; `singularity pull`
   CellBender `.sif`; set `env/paths.sh`.
3. **Build `manifest.tsv`** from staged raw matrices.
4. **`bash slurm/submit_all.sh`** — launches the dependency-chained pipeline (open arms).
5. Monitor: `squeue -u $USER`, `sacct -j <jid> --format=JobID,State,Elapsed,MaxRSS`, logs in `$PROJ/logs`.
6. After QC gate + DE, **review the two gates** (decontamination dual-count concordance; procurement
   robustness) before trusting any DEG — exactly as in v2.1.
7. Fold in controlled arms only if `access_status.md` flips to `raw-reprocess`; re-run Steps 2→ for
   those samples via the array.

---

*Generated as a BigPurple adaptation of Plan v2.1. Scientific design, contrasts, gates, and gene
expectations are carried forward unchanged; only the compute/orchestration layer is new.*
