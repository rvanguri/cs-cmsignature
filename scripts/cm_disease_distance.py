#!/usr/bin/env python3
"""Do CS cardiomyocytes cluster apart from cardiomyocytes of other cardiomyopathies?

Cell-intrinsic specificity figure for the CM-core story. Reads the per-patient cardiomyocyte
pseudobulk (cm_pseudobulk_counts.tsv genes x patients + cm_pseudobulk_meta.tsv), CPM-log normalizes,
and asks whether the CS cardiomyocyte profile is a distinct state:

  - PCA of patients coloured by DISEASE, and a parallel panel coloured by STUDY (CS = one cohort, so
    the study panel keeps the confound visible instead of hiding it).
  - disease-averaged correlation heatmap.
  - centroid-distance table: is CS farther from each mimic than the mimics are from each other?
  - silhouette of CS vs non-CS.

Runs on all genes (top-variance) by default; pass --signature to restrict to the validated CS-CM set
and show separation on the signature specifically.
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import log, require, ensure_dir  # noqa: E402
LOG = log("cm_distance")

DISEASE_COLORS = {"CS": "#d62728", "DCM": "#1f77b4", "ARVC": "#2ca02c", "HCM": "#9467bd",
                  "NF": "#7f7f7f", "ICM": "#ff7f0e", "NCC": "#8c564b"}
STUDY_MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pb_counts", required=True)
    ap.add_argument("--pb_meta", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--signature", default=None, help="optional 1-col gene list to restrict to")
    ap.add_argument("--n_hvg", type=int, default=2000)
    ap.add_argument("--min_mean_cpm", type=float, default=1.0)
    args = ap.parse_args()
    require(args.pb_counts, "pb counts"); require(args.pb_meta, "pb meta"); ensure_dir(args.out)

    from sklearn.decomposition import PCA
    from sklearn.metrics import silhouette_score
    from scipy.spatial.distance import pdist, squareform
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    counts = pd.read_csv(args.pb_counts, sep="\t", index_col=0)
    meta = pd.read_csv(args.pb_meta, sep="\t", index_col=0)
    meta = meta.reindex(counts.columns)
    disease = meta["disease"].astype(str)
    study = meta.get("study", pd.Series("NA", index=meta.index)).astype(str)
    LOG.info("patients=%d genes=%d diseases=%s", counts.shape[1], counts.shape[0],
             dict(disease.value_counts()))

    # CPM + log1p
    cpm = counts.divide(counts.sum(0), axis=1) * 1e6
    logcpm = np.log1p(cpm)
    logcpm = logcpm[(cpm.mean(1) >= args.min_mean_cpm)]           # drop low-expression genes

    if args.signature and os.path.exists(args.signature):
        sig = pd.read_csv(args.signature, sep="\t", header=None)[0].astype(str).tolist()
        genes = [g for g in sig if g in logcpm.index]
        LOG.info("restricting to %d/%d signature genes", len(genes), len(sig))
        X = logcpm.loc[genes]
        tag = "signature"
    else:
        v = logcpm.var(1).sort_values(ascending=False)
        X = logcpm.loc[v.index[:args.n_hvg]]
        tag = f"hvg{args.n_hvg}"

    # z-score genes, patients x genes
    Z = ((X.T - X.T.mean()) / (X.T.std() + 1e-9)).fillna(0.0)
    pca = PCA(n_components=min(10, Z.shape[0] - 1, Z.shape[1])).fit(Z)
    pcs = pd.DataFrame(pca.transform(Z), index=Z.index,
                       columns=[f"PC{i+1}" for i in range(pca.n_components_)])
    ev = pca.explained_variance_ratio_ * 100

    # ---- PCA scatter: by disease, and by study ----
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    studies = sorted(study.unique())
    smark = {s: STUDY_MARKERS[i % len(STUDY_MARKERS)] for i, s in enumerate(studies)}
    for ax, colour_by, title in [(axes[0], "disease", "coloured by disease"),
                                 (axes[1], "study", "coloured by study")]:
        for p in pcs.index:
            d, s = disease[p], study[p]
            if colour_by == "disease":
                c = DISEASE_COLORS.get(d, "#333333"); mk = smark[s]
            else:
                pal = plt.cm.tab10(np.linspace(0, 1, len(studies)))
                c = pal[studies.index(s)]; mk = "o"
            ax.scatter(pcs.loc[p, "PC1"], pcs.loc[p, "PC2"], c=[c], marker=mk, s=55,
                       edgecolors="k", linewidths=0.3, alpha=0.85)
        ax.set_xlabel(f"PC1 ({ev[0]:.0f}%)"); ax.set_ylabel(f"PC2 ({ev[1]:.0f}%)")
        ax.set_title(f"CM pseudobulk PCA ({tag}) — {title}")
    # legends
    from matplotlib.lines import Line2D
    dl = [Line2D([0], [0], marker="o", color="w", markerfacecolor=DISEASE_COLORS.get(d, "#333"),
                 markeredgecolor="k", markersize=8, label=d) for d in disease.unique()]
    axes[0].legend(handles=dl, fontsize=8, title="disease", loc="best")
    fig.tight_layout(); fig.savefig(os.path.join(args.out, f"cm_pca_{tag}.png"), dpi=200,
                                    bbox_inches="tight"); plt.close(fig)

    # ---- disease-averaged correlation heatmap ----
    prof = X.T.groupby(disease.values).mean().T            # genes x disease
    corr = prof.corr(method="spearman")
    fig, ax = plt.subplots(figsize=(1.1 * len(corr) + 2, 1.0 * len(corr) + 1))
    im = ax.imshow(corr.values, vmin=corr.values.min(), vmax=1, cmap="magma")
    ax.set_xticks(range(len(corr))); ax.set_xticklabels(corr.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(corr))); ax.set_yticklabels(corr.index)
    for i in range(len(corr)):
        for j in range(len(corr)):
            ax.text(j, i, f"{corr.values[i, j]:.2f}", ha="center", va="center",
                    color="w" if corr.values[i, j] < 0.9 else "k", fontsize=8)
    fig.colorbar(im, ax=ax, label="Spearman r (disease-mean CM profile)")
    fig.tight_layout(); fig.savefig(os.path.join(args.out, f"cm_disease_corr_{tag}.png"), dpi=200,
                                    bbox_inches="tight"); plt.close(fig)

    # ---- centroid distances (in PCA space) ----
    cent = pcs.groupby(disease.values).mean()
    D = pd.DataFrame(squareform(pdist(cent.values)), index=cent.index, columns=cent.index)
    D.to_csv(os.path.join(args.out, f"cm_disease_distances_{tag}.tsv"), sep="\t")
    if "CS" in D.index:
        cs_d = D.loc["CS"].drop("CS").sort_values()
        others = D.drop(index="CS", columns="CS")
        mean_between_mimics = others.values[np.triu_indices(len(others), 1)].mean()
        LOG.info("CS centroid distance to each mimic:\n%s", cs_d.round(2).to_string())
        LOG.info("mean CS->mimic distance = %.2f | mean mimic<->mimic distance = %.2f",
                 cs_d.mean(), mean_between_mimics)

    # ---- silhouette CS vs non-CS ----
    labels = (disease == "CS").astype(int).values
    if len(set(labels)) == 2:
        sil = silhouette_score(pcs.values, labels)
        LOG.info("silhouette (CS vs non-CS) in PCA space = %.3f", sil)
    pcs.assign(disease=disease.values, study=study.values).to_csv(
        os.path.join(args.out, f"cm_pca_coords_{tag}.tsv"), sep="\t")
    LOG.info("done -> %s", args.out)


if __name__ == "__main__":
    main()
