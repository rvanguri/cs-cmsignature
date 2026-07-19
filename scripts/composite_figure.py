#!/usr/bin/env python3
"""Assemble the 6-panel Research Letter figure from local analysis tables (draft for Illustrator polish).

A  PCA of per-patient cardiomyocyte pseudobulk (CS vs mimics)
B  four-gene signature logFC across CS-vs-mimic contrasts + Foong spatial replication (heatmap)
C  content-normalized recovery in CM-dominant Visium spots, CS vs normal (grouped bars, fold labels)
D  signature scores by CS tissue zone (CM-intrinsic vs curated granuloma-inflammation reference)
E  CM-intrinsic vs distance-to-lesion within preserved spots (binned mean + mixed-model annotation)
F  representative in-situ H&E-paired CM-intrinsic map (embedded)
"""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from sklearn.decomposition import PCA

BASE = os.environ.get("PROJECT_ROOT", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(BASE, "figures", "Figure1_CS_cardiomyocyte.png")
os.makedirs(os.path.dirname(OUT), exist_ok=True)
GENES = ["GJB7", "TNNI3K", "MLIP", "PANK1"]
DCOL = {"CS":"#d62728","DCM":"#1f77b4","ARVC":"#2ca02c","HCM":"#9467bd","NF":"#7f7f7f","ICM":"#ff7f0e","NCC":"#8c564b"}
plt.rcParams.update({"font.size":8,"font.family":"DejaVu Sans","axes.titlesize":9})

def panel_label(ax, L):
    ax.text(-0.12, 1.06, L, transform=ax.transAxes, fontsize=13, fontweight="bold", va="top", ha="right")

fig = plt.figure(figsize=(11, 9.4))
# 3 rows: [A B C] / [D E] / [F full-width banner] so the F in-situ strip gets the whole width
gs = GridSpec(3, 3, figure=fig, height_ratios=[1.0, 1.0, 0.95], hspace=0.5, wspace=0.38)

# ---- A: PCA ----
axA = fig.add_subplot(gs[0,0])
counts = pd.read_csv(f"{BASE}/results/de/cm_pseudobulk_counts.tsv", sep="\t", index_col=0)
meta = pd.read_csv(f"{BASE}/results/de/cm_pseudobulk_meta.tsv", sep="\t", index_col=0).reindex(counts.columns)
cpm = counts.divide(counts.sum(0), axis=1)*1e6; logcpm = np.log1p(cpm)
logcpm = logcpm[cpm.mean(1) >= 1.0]
hv = logcpm.var(1).sort_values(ascending=False).index[:2000]
Z = ((logcpm.loc[hv].T - logcpm.loc[hv].T.mean())/(logcpm.loc[hv].T.std()+1e-9)).fillna(0)
pcs = PCA(5).fit(Z); XY = pcs.transform(Z); ev = pcs.explained_variance_ratio_*100
dis = meta["disease"].astype(str).values
for d in pd.unique(dis):
    m = dis==d
    axA.scatter(XY[m,0], XY[m,1], s=22, c=DCOL.get(d,"#333"), edgecolors="k", linewidths=0.3,
                alpha=0.85, label=d)
axA.set_xlabel(f"PC1 ({ev[0]:.0f}%)"); axA.set_ylabel(f"PC2 ({ev[1]:.0f}%)")
axA.set_title("Cardiomyocyte pseudobulk"); axA.legend(fontsize=6, ncol=2, frameon=False, loc="best")
panel_label(axA, "A")

# ---- B: signature heatmap ----
axB = fig.add_subplot(gs[0,1])
val = pd.read_csv(f"{BASE}/results/de/de_cm/CM_foong_validation.tsv", sep="\t", index_col=0)
cols = ["lfc_DCM","lfc_ARVC","lfc_HCM","lfc_NF","foong_log2FC"]
H = val.loc[GENES, cols].astype(float)
H.columns = ["vs DCM","vs ARVC","vs HCM","vs NF","Visium\n(vs HCM)"]
vmax = np.nanpercentile(np.abs(H.values), 98)
im = axB.imshow(H.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
axB.set_xticks(range(len(H.columns))); axB.set_xticklabels(H.columns, rotation=45, ha="right", fontsize=7)
axB.set_yticks(range(len(GENES))); axB.set_yticklabels(GENES, style="italic")
for i in range(len(GENES)):
    for j in range(len(H.columns)):
        axB.text(j, i, f"{H.values[i,j]:+.1f}", ha="center", va="center", fontsize=6.5,
                 color="white" if abs(H.values[i,j])>vmax*0.6 else "black")
axB.set_title("Signature log$_2$FC (CS up)")
fig.colorbar(im, ax=axB, fraction=0.046, pad=0.04, label="log$_2$FC")
panel_label(axB, "B")

# ---- C: content-normalized recovery ----
axC = fig.add_subplot(gs[0,2])
ps = pd.read_csv(f"{BASE}/results/figures/cm_content_normalized/cm_content_perspot.tsv.gz", sep="\t")
mean = ps.groupby(["gene","condition"])["expr"].mean().unstack()
mean = mean.reindex(GENES)
x = np.arange(len(GENES)); w = 0.38
axC.bar(x-w/2, mean["CS"], w, color="#d62728", label="CS")
axC.bar(x+w/2, mean["control"], w, color="#9aa0a6", label="normal")
ymax = float(mean[["CS","control"]].max().max())
axC.set_ylim(top=ymax*1.12)                       # headroom so fold labels clear the top frame
for i,g in enumerate(GENES):
    cs, ct = mean.loc[g,"CS"], mean.loc[g,"control"]
    fold = "n.d." if ct<=1e-6 else f"{cs/ct:.1f}×"
    axC.text(i, cs + ymax*0.015, fold, ha="center", va="bottom", fontsize=6.5)
axC.set_xticks(x); axC.set_xticklabels(GENES, style="italic", rotation=30, ha="right")
axC.set_ylabel("mean expr (CM-dominant spots)"); axC.set_title("Visium CS vs normal")
axC.legend(fontsize=6, frameon=False)
panel_label(axC, "C")

# ---- D: zones ----
# Center D and E across the full middle row (margin columns on each side) so they
# don't sit left-justified with an empty third cell.
gsDE = GridSpecFromSubplotSpec(1, 4, subplot_spec=gs[1, :],
                               width_ratios=[0.5, 1.0, 1.0, 0.5], wspace=0.42)
axD = fig.add_subplot(gsDE[0, 1])
pz = pd.read_csv(f"{BASE}/results/figures/foong_regional/foong_regional_perspot.tsv.gz", sep="\t")
zorder = [z for z in ["preserved","granulomatous","fibrotic"] if z in pz.zone.unique()]
mI = [pz.loc[pz.zone==z,"CM_intrinsic"].mean() for z in zorder]
sI = [pz.loc[pz.zone==z,"CM_intrinsic"].sem() for z in zorder]
mF = [pz.loc[pz.zone==z,"inflammatory"].mean() for z in zorder]
sF = [pz.loc[pz.zone==z,"inflammatory"].sem() for z in zorder]
x = np.arange(len(zorder)); w=0.38
axD.bar(x-w/2, mI, w, yerr=sI, capsize=2, color="#1b9e77", label="CM-intrinsic")
axD.bar(x+w/2, mF, w, yerr=sF, capsize=2, color="#d95f02", label="granuloma inflammation")
axD.axhline(0, color="k", lw=0.5)
axD.set_xticks(x); axD.set_xticklabels(zorder, rotation=30, ha="right"); axD.set_ylabel("signature score")
axD.set_title("Visium signature by CS tissue zone"); axD.legend(fontsize=6, frameon=False)
panel_label(axD, "D")

# ---- E: distance to lesion ----
axE = fig.add_subplot(gsDE[0, 2])
_dpath = f"{BASE}/results/figures/foong_regional/foong_distance_to_lesion.tsv.gz"
dd = None
try:
    if os.path.getsize(_dpath) > 0:
        dd = pd.read_csv(_dpath, sep="\t")
except Exception:
    dd = None
if dd is None or dd.empty:
    axE.text(0.5, 0.5, "[panel E: re-stage\nfoong_distance_to_lesion.tsv.gz]",
             ha="center", va="center", fontsize=8, color="#b00")
else:
    dd["dz"] = dd.groupby("sample")["dist"].transform(lambda x:(x-x.mean())/(x.std()+1e-9))
    dd = dd.replace([np.inf,-np.inf],np.nan).dropna(subset=["dz","CM_intrinsic"])
    dd["bin"] = pd.qcut(dd["dz"].rank(method="first"), 10, labels=False)
    cv = dd.groupby("bin").agg(d=("dz","mean"), s=("CM_intrinsic","mean"),
                               e=("CM_intrinsic",lambda x:x.std()/np.sqrt(len(x)))).reset_index()
    axE.errorbar(cv["d"], cv["s"], yerr=cv["e"], marker="o", ms=4, capsize=2, color="#333",
                 ls="none", label="distance-decile mean ± SEM (10 bins)")
    b = np.polyfit(dd["dz"], dd["CM_intrinsic"], 1)
    xs = np.linspace(dd["dz"].min(), dd["dz"].max(), 40)
    axE.plot(xs, b[0]*xs+b[1], "r-", lw=1.5, label="spot-level linear fit")
    axE.legend(fontsize=5.5, frameon=False, loc="upper left")
axE.set_xlabel("distance to nearest lesion (within-sample z)"); axE.set_ylabel("CM-intrinsic score")
axE.set_title("positive slope β=+0.044/SD (P<0.001),\npositive in 8/8 patients", fontsize=7.5)
panel_label(axE, "E")

# ---- F: representative spatial H&E-paired map (full-width bottom banner) ----
axF = fig.add_subplot(gs[2, :]); axF.axis("off")
# 3-panel strip: H&E | granuloma (immune) score | CM-intrinsic score, for the representative CS sample.
fpath = f"{BASE}/results/figures/panelF_candidates/panelF_CS_104.png"
if os.path.exists(fpath):
    raw = mpimg.imread(fpath)
    # Crop off the baked-in per-subpanel sample titles (top band) and white margins so the three
    # in-situ maps fill the banner. Trim a fixed top title band, then tight-crop to non-white content.
    h, w = raw.shape[:2]
    band = int(round(h * 0.146))               # top title band (~"CS_104 ..." captions)
    body = raw[band:, :, :]
    nonwhite = (body[:, :, :3].min(axis=2) < 0.96)
    rows = np.where(nonwhite.mean(axis=1) > 0.01)[0]
    cols = np.where(nonwhite.mean(axis=0) > 0.01)[0]
    if rows.size and cols.size:
        body = body[rows.min():rows.max()+1, cols.min():cols.max()+1, :]
    axF.imshow(body, interpolation="lanczos")
    axF.set_title("Visium: H&E | immune-cell density (granuloma) | CM-intrinsic score", fontsize=8)
else:
    axF.text(0.5,0.5,"[run foong_panelF_candidates.py]",ha="center",va="center")
panel_label(axF, "F")

fig.savefig(OUT, dpi=600, bbox_inches="tight")
fig.savefig(OUT.replace(".png",".pdf"), bbox_inches="tight")
print("wrote", OUT)
