#!/usr/bin/env python3
"""Two-part figure for the JACC Advances Brief Report (max 1 simple figure, <=2 parts).

Part A  four-gene CS-cardiomyocyte signature log2FC across CS-vs-mimic pseudobulk contrasts
        (DCM/ARVC/HCM/NF) + independent Foong Visium spatial replication column (heatmap).
Part B  representative in-situ Visium strip for CS_104: H&E | immune-cell density (granuloma)
        | CM-intrinsic signature score.

Layout: Part A (compact heatmap, centered) stacked above Part B (full-width in-situ banner).
Both parts carry spatial data: A includes the Visium replication column; B is spatial in situ.
"""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

BASE = os.environ.get("PROJECT_ROOT", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(BASE, "figures", "Figure1_JACC_briefreport.png")
os.makedirs(os.path.dirname(OUT), exist_ok=True)
GENES = ["GJB7", "TNNI3K", "MLIP", "PANK1"]
plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans", "axes.titlesize": 9})

def part_label(ax, L):
    ax.text(-0.12, 1.06, L, transform=ax.transAxes, fontsize=14, fontweight="bold",
            va="top", ha="right")

fig = plt.figure(figsize=(7.2, 7.0))
gs = GridSpec(2, 1, figure=fig, height_ratios=[1.0, 1.05], hspace=0.32)

# ---- A: signature heatmap (centered in the top row via margin columns) ----
gsA = GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[0, 0],
                              width_ratios=[0.55, 1.0, 0.55], wspace=0.0)
axA = fig.add_subplot(gsA[0, 1])
val = pd.read_csv(f"{BASE}/results/de/de_cm/CM_foong_validation.tsv", sep="\t", index_col=0)
cols = ["lfc_DCM", "lfc_ARVC", "lfc_HCM", "lfc_NF", "foong_log2FC"]
H = val.loc[GENES, cols].astype(float)
H.columns = ["vs DCM", "vs ARVC", "vs HCM", "vs NF", "Visium\n(vs HCM)"]
vmax = np.nanpercentile(np.abs(H.values), 98)
im = axA.imshow(H.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
axA.set_xticks(range(len(H.columns)))
axA.set_xticklabels(H.columns, rotation=45, ha="right", fontsize=7.5)
axA.set_yticks(range(len(GENES)))
axA.set_yticklabels(GENES, style="italic", fontsize=9)
for i in range(len(GENES)):
    for j in range(len(H.columns)):
        axA.text(j, i, f"{H.values[i, j]:+.1f}", ha="center", va="center", fontsize=7.5,
                 color="white" if abs(H.values[i, j]) > vmax * 0.6 else "black")
axA.set_title("CS-cardiomyocyte signature log$_2$FC (CS up)", fontsize=9)
fig.colorbar(im, ax=axA, fraction=0.046, pad=0.04, label="log$_2$FC")
part_label(axA, "A")

# ---- B: representative in-situ Visium strip (full-width banner) ----
axB = fig.add_subplot(gs[1, 0]); axB.axis("off")
fpath = f"{BASE}/results/figures/panelF_candidates/panelF_CS_104.png"
if os.path.exists(fpath):
    raw = mpimg.imread(fpath)
    h, w = raw.shape[:2]
    band = int(round(h * 0.146))               # trim baked per-subpanel sample captions
    body = raw[band:, :, :]
    nonwhite = (body[:, :, :3].min(axis=2) < 0.96)
    rows = np.where(nonwhite.mean(axis=1) > 0.01)[0]
    cols_ = np.where(nonwhite.mean(axis=0) > 0.01)[0]
    if rows.size and cols_.size:
        body = body[rows.min():rows.max() + 1, cols_.min():cols_.max() + 1, :]
    axB.imshow(body, interpolation="lanczos")
    axB.set_title("Visium in situ: H&E | immune-cell density (granuloma) | CM-intrinsic score",
                  fontsize=8.5)
else:
    axB.text(0.5, 0.5, "[run foong_panelF_candidates.py]", ha="center", va="center")
part_label(axB, "B")

fig.savefig(OUT, dpi=600, bbox_inches="tight")
fig.savefig(OUT.replace(".png", ".pdf"), bbox_inches="tight")
print("wrote", OUT)
