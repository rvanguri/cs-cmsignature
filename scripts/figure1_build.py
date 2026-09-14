#!/usr/bin/env python3
"""Build the manuscript figure from the committed panel value tables.

Panel a: the 134 core genes testable in both arms -- spatial log2FC against pooled
single-nucleus log2FC. Quadrant shading marks sign concordance. No regression line is
drawn: the two axes are different assays on different cohorts and are not calibrated to
each other, so a fitted slope would invite a magnitude reading the data do not support.
Genes failing locus-level QC are drawn as open markers and still plotted, because they
are part of the tested set; they are excluded from no statistic shown here.

Panel b: the 18 genes concordant in direction and significant in both arms, with the
two published cardiomyocyte comparisons against non-failing donors alongside. The
published tables list only genes called differential, so a gene absent from one is drawn
as an open square -- not as a null effect. Colour scales differ between columns because
the assays have different dynamic ranges; that is stated in the colorbar label.

Inputs   results/tables/figure1_panelA_values.tsv
         results/tables/figure1_panelB_values.tsv
Outputs  figures/Figure1.png, figures/Figure1.pdf
"""
import argparse

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from adjustText import adjust_text
from matplotlib.colors import SymLogNorm
from matplotlib.lines import Line2D
from scipy import stats

sys_kernel = None
try:  # figure-style kernel helpers when available
    apply_figure_style  # type: ignore[name-defined]
except NameError:
    from importlib import import_module
    try:
        sys_kernel = import_module("figure_style_kernel")
    except ModuleNotFoundError:
        sys_kernel = None

TAB = "results/tables"
BLUE, RED, GREY = "#2166ac", "#b2182b", "#9a9a9a"
CMAP = "RdBu_r"
# Panel b colour: one shared scale for all four columns, symmetric about zero, with a
# symmetric-log transform so the single-nucleus effects (|log2FC| < 1, the biological
# scale of a pseudobulk contrast) stay legible on the same key as the spatial effects
# (|log2FC| up to 6, inflated by probe-level quantification). Ticks are placed at the
# breakpoints of that transform, not at even spacing.
TICKS = [-6, -3, -1, 0, 1, 3, 6]
LINTHRESH = 1.0
COLS = [("logFC", "adj.P.Val", "spatial"), ("lfc_pool", "fdr_pool", "snRNA"),
        ("lfc_DCMNF", "fdr_DCMNF", "DCM\nvs NF"), ("lfc_HCMNF", "fdr_HCMNF", "HCM\nvs NF")]


def fdr_area(fdr, lo=18, hi=190):
    """Marker area scaled by -log10 FDR, saturating at 1e-4."""
    z = np.clip(-np.log10(np.clip(fdr, 1e-300, 1)), 0, 4) / 4
    return lo + (hi - lo) * z


def panel_a(ax, pA):
    rho, p = stats.spearmanr(pA.spatial_logFC, pA.lfc_pool)
    # One gene (HSPB2, pooled log2FC +7.1) sits an order of magnitude above the rest of
    # the single-nucleus effects; a linear y-axis would leave 80% of the range empty.
    # A symmetric-log y keeps every point plotted and readable, with ticks at the
    # breakpoints so the compression is visible rather than implied.
    ax.set_yscale("symlog", linthresh=1.0, linscale=0.75)
    ax.set_yticks([-2, -1, 0, 1, 2, 5])
    ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.yaxis.set_minor_locator(mpl.ticker.NullLocator())
    xlo, xhi = pA.spatial_logFC.min() - 0.6, pA.spatial_logFC.max() + 0.9
    ylo, yhi = -2.6, 8.5
    # concordant quadrants: both negative, or both positive
    ax.add_patch(plt.Rectangle((xlo, ylo), -xlo, -ylo, fc="#dfe6ee", ec="none", zorder=0))
    ax.add_patch(plt.Rectangle((0, 0), xhi, yhi, fc="#dfe6ee", ec="none", zorder=0))
    ax.axhline(0, color=GREY, lw=0.6, zorder=1)
    ax.axvline(0, color=GREY, lw=0.6, zorder=1)

    for direction, colour in (("down", BLUE), ("up", RED)):
        for qc, kw in ((True, dict(facecolors=colour, edgecolors="none")),
                       (False, dict(facecolors="none", edgecolors=colour, linewidths=0.9))):
            s = pA[(pA.direction == direction) & (pA.qc_pass.astype(bool) == qc)]
            ax.scatter(s.spatial_logFC, s.lfc_pool, s=26, zorder=3, **kw)

    ax.set_xlim(xlo, xhi)
    ax.set_ylim(ylo, yhi)

    lab = pA[pA.label_gene.astype(bool)]
    texts = [ax.text(r.spatial_logFC, r.lfc_pool, f"$\\it{{{r.gene}}}$", fontsize=6,
                     color=BLUE if r.direction == "down" else RED, zorder=4)
             for _, r in lab.iterrows()]
    adjust_text(texts, x=lab.spatial_logFC.values, y=lab.lfc_pool.values, ax=ax,
                expand=(1.35, 1.6), force_text=(0.4, 0.6), force_static=(0.15, 0.25),
                min_arrow_len=3, avoid_self=False, time_lim=8, ensure_inside_axes=True,
                arrowprops=dict(arrowstyle="-", lw=0.5, color=GREY, shrinkA=1, shrinkB=2))

    ax.text(0.02, 0.97, f"Spearman $\\rho$ = {rho:+.2f} (n = {len(pA)}, p = {p:.2f})\n"
                        "no line fitted: axes are not calibrated to each other",
            transform=ax.transAxes, va="top", ha="left", fontsize=6)
    ax.text(0.97, 0.80, "concordant", transform=ax.transAxes, ha="right",
            fontsize=6, color=GREY)
    ax.text(0.03, 0.06, "concordant", transform=ax.transAxes, ha="left",
            fontsize=6, color=GREY)
    ax.set_xlabel("spatial log$_2$FC, CS vs comparator cardiomyopathies")
    ax.set_ylabel("single-nucleus log$_2$FC,\nCS vs DCM + ARVC (pooled)")
    ax.set_title("Direction replicates across platforms; magnitude does not", fontsize=8)
    ax.legend(handles=[
        Line2D([], [], marker="o", ls="", mfc=BLUE, mec="none", ms=5, label="down in CS (spatial)"),
        Line2D([], [], marker="o", ls="", mfc=RED, mec="none", ms=5, label="up in CS (spatial)"),
        Line2D([], [], marker="o", ls="", mfc="none", mec=GREY, ms=5, label="fails locus QC")],
        loc="lower right", frameon=False, fontsize=6, handletextpad=0.4, borderpad=0.2)
    return rho, p


def panel_b(ax, pB, cax_col, cax_size):
    pB = pB.sort_values(["direction", "spatial_logFC"], ascending=[True, True]).reset_index(drop=True)
    n = len(pB)
    y = np.arange(n)[::-1]
    norm = SymLogNorm(linthresh=LINTHRESH, vmin=-6, vmax=6, base=10)
    cmap = plt.get_cmap(CMAP)

    for j, (lfc_c, fdr_c, label) in enumerate(COLS):
        v, f = pB[lfc_c].values, pB[fdr_c].values
        drawn = ~pd.isna(v)
        ax.scatter(np.full(drawn.sum(), j), y[drawn], s=fdr_area(f[drawn]),
                   c=cmap(norm(v[drawn])), edgecolors="#4d4d4d", linewidths=0.4, zorder=3)
        ax.scatter(np.full((~drawn).sum(), j), y[~drawn], s=44, marker="s",
                   facecolors="none", edgecolors=GREY, linewidths=0.7, zorder=3)

    # separate the down- and up-in-CS blocks
    split = int((pB.direction == "down").sum())
    ax.axhline(y[split] + 0.5, color="#4d4d4d", lw=0.7)
    for edge in (y.max() + 0.5, y.min() - 0.5):
        ax.axhline(edge, color="#4d4d4d", lw=0.7)

    marks = {"locus_quant_outlier": "$^\\dagger$", "no_annotation": "$^{*}$"}
    ax.set_yticks(y)
    ax.set_yticklabels([f"$\\it{{{g}}}${marks.get(m, '')}"
                        for g, m in zip(pB.gene, pB.qc_mode)])
    ax.set_xticks(range(len(COLS)))
    ax.set_xticklabels([c[2] for c in COLS])
    ax.set_xlim(-0.6, len(COLS) - 0.4)
    ax.set_ylim(y.min() - 0.9, y.max() + 1.0)
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(False)
    ax.tick_params(length=0)
    ax.set_title("Reduced in cardiomyopathy, reduced further in CS", fontsize=8)

    # group headers over the two assay families
    for x0, x1, text in ((0, 1, "this study"), (2, 3, "published (Chaffin et al.)")):
        ax.plot([x0 - 0.35, x1 + 0.35], [y.max() + 0.75] * 2, color="#4d4d4d",
                lw=0.7, clip_on=False)
        ax.text((x0 + x1) / 2, y.max() + 0.95, text, ha="center", va="bottom", fontsize=6)

    cb = mpl.colorbar.Colorbar(cax_col, cmap=cmap, norm=norm,
                               orientation="horizontal", ticks=TICKS)
    cb.set_label("log$_2$FC", fontsize=6, labelpad=1)
    cb.ax.tick_params(labelsize=6, length=2, pad=1)
    cb.ax.xaxis.set_minor_locator(mpl.ticker.NullLocator())
    cb.ax.xaxis.set_major_formatter(mpl.ticker.FixedFormatter([f"{t:g}" for t in TICKS]))
    cb.outline.set_visible(False)

    cax_size.set_axis_off()
    cax_size.set_xlim(0, 1)
    cax_size.set_ylim(0, 1)
    for k, (fdr, lbl) in enumerate([(0.05, "0.05"), (0.01, "0.01"), (1e-4, "$\\leq 10^{-4}$")]):
        cax_size.scatter([0.16], [0.62 - 0.27 * k], s=fdr_area(fdr), c="#d9d9d9",
                         edgecolors="#4d4d4d", linewidths=0.4, clip_on=False)
        cax_size.text(0.42, 0.62 - 0.27 * k, lbl, va="center", fontsize=6)
    cax_size.text(0.0, 0.95, "FDR", ha="left", va="center", fontsize=6)
    return pB


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="figures/Figure1")
    args = ap.parse_args()
    if sys_kernel is None:
        mpl.rcParams.update({"font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 6,
                             "ytick.labelsize": 6, "axes.spines.top": False,
                             "axes.spines.right": False, "font.family": "sans-serif",
                             "pdf.fonttype": 42})

    pA = pd.read_csv(f"{TAB}/figure1_panelA_values.tsv", sep="\t")
    pB = pd.read_csv(f"{TAB}/figure1_panelB_values.tsv", sep="\t")

    fig = plt.figure(figsize=(9.2, 5.6))
    axA = fig.add_axes([0.065, 0.115, 0.44, 0.80])
    axB = fig.add_axes([0.635, 0.215, 0.20, 0.70])
    cax_col = fig.add_axes([0.615, 0.085, 0.17, 0.020])
    cax_size = fig.add_axes([0.865, 0.055, 0.11, 0.13])

    rho, p = panel_a(axA, pA)
    pBs = panel_b(axB, pB, cax_col, cax_size)

    fig.text(0.565, 0.020, "$^\\dagger$locus-quantification outlier   $^{*}$no annotation   "
                           "open square: gene not reported in that table", fontsize=6)
    for ax, x in ((axA, 0.012), (axB, 0.545)):
        fig.text(x, 0.955, "a" if ax is axA else "b", fontsize=11,
                 fontweight="bold", va="top", ha="left")

    fig.savefig(f"{args.out}.png", dpi=300)
    fig.savefig(f"{args.out}.pdf")
    n_open = int((~pA.qc_pass.astype(bool)).sum())
    print(f"[panel a] n={len(pA)} | rho={rho:+.3f} p={p:.3f} | open markers {n_open}")
    print(f"[panel b] n={len(pBs)} | "
          f"{int((pBs.direction == 'down').sum())} down, {int((pBs.direction == 'up').sum())} up")
    for lfc_c, fdr_c, label in COLS[2:]:
        print(f"[panel b] {label.replace(chr(10), ' ')}: "
              f"{int(pBs[lfc_c].notna().sum())} of {len(pBs)} reported")


if __name__ == "__main__":
    main()
