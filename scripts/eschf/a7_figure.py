#!/usr/bin/env python3
"""ESC HF A7: four-panel Figure 1, graphical abstract, legend (revision v3: HANDOFF_figure_v3_ESCHF.md).

Run from repo root after the ESCHF_PRIMARY=v2 cascade, a2_perpatient.py and a9_snrna_screen.py. All panels come
from the primary fit over all cardiomyocyte-dominant spots (results/eschf/primary_v2/*). v3 changes: per-patient
(old D) and Hamed (old F) panels removed, 2 x 2 layout; Panel B relabelled ("Score of genes lower in CS"), slope
annotations above the axes; Panels C and D apply snRNA screen 3a (|pooled snRNA log2FC| > 3; a9_snrna_screen.tsv);
Panel D = concordant genes from all spatial hits (not the section-stable core), with corrected locus-QC symbols.
v4 changes (user instructions): Panel C screened genes (3a failures) drawn as open neutral-grey triangles with no
sign category; Panel B legend pathway sentence replaced (asserted against a14_pathway_hits_removed.tsv); Panel D =
the concordant genes of the section-stable core after screen 3a (in_concordant_core_v2_screened, 17 genes;
NDUFC2 removed, v5). Panel C legend states the 3b source-absence pattern of the screened genes (asserted).
"Granuloma" in earlier versions = immune-dominant spot (column and file names unchanged).
Every plotted number is written to results/eschf/a7_panel_values.tsv for manuscript_results.py. H&E images
(CS and comparator) are the downscaled JPEGs in data/eschf/img, with scale factors in data/eschf/img/image_meta.json.
"""
import glob, json
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from adjustText import adjust_text
from matplotlib.collections import EllipseCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from PIL import Image
from scipy import stats

OUT, TOP, FIG = "results/eschf/primary_v2", "results/eschf", "figures"
YLIM_C = 3.0    # Panel C y-axis half-range = SCREEN; points beyond are drawn at the edge as triangles
SCREEN = 3.0    # snRNA screen 3a: |pooled snRNA log2FC| > 3 is unreliable (pre-specified)
CS_SEC, COMP_SEC = "CS_108-1", "ARVC_1"
mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "font.size": 7, "axes.titlesize": 7, "axes.labelsize": 7, "xtick.labelsize": 6,
                     "ytick.labelsize": 6, "legend.fontsize": 6, "axes.spines.top": False,
                     "axes.spines.right": False, "pdf.fonttype": 42, "svg.fonttype": "none",
                     "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6})
# Okabe-Ito
C_CS, C_CMP = "#D55E00", "#0072B2"
C_SEL, C_NEAR, C_IMM, C_OTH = "#009E73", "#9AD9C2", "#E69F00", "#C8C8C8"
C_BG = "#BDBDBD"
C_SCR = "0.35"   # screened genes (Panel C): open triangles, neutral grey
C_LO, C_HI = "#0072B2", "#E69F00"
VALS = []


def val(panel, name, value):
    VALS.append(dict(panel=panel, quantity=name, value=value))


def it(g):
    return r"$\it{" + g + "}$"


def letter(ax, s, dx=-0.02, dy=1.02):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=9, fontweight="bold", ha="right", va="bottom")


# ------------------------------------------------------------------ data
sp = pd.read_csv("data/eschf/primary_v2/eschf_perspot_T2000.tsv.gz", sep="\t")
imeta = json.load(open("data/eschf/img/image_meta.json"))
bins = pd.read_csv(f"{OUT}/a2_distance_bins.tsv", sep="\t")
bins = bins[bins.n_spots > 0].reset_index(drop=True)   # populated bins only
slp = pd.read_csv(f"{TOP}/a2_perpatient_slopes.tsv", sep="\t").set_index("outcome")
bpp = pd.read_csv(f"{OUT}/a2_distance_bins_perpatient.tsv", sep="\t")
cref = pd.read_csv(f"{OUT}/a2_comparator_reference.tsv", sep="\t").set_index("score")
H = pd.read_csv(f"{OUT}/q1_hits_annotated.tsv", sep="\t").set_index("gene")
SN = pd.read_csv(f"{OUT}/snrna_pooled_full.tsv", sep="\t").set_index("gene")
NF = pd.read_csv("data/snrna/snrna_CS_vs_NF.tsv", sep="\t").set_index("gene")
pmap = pd.read_csv(f"{OUT}/section_patient_map.tsv", sep="\t")
qc = pd.read_csv(f"{OUT}/a3_qc_comparison.tsv", sep="\t")
SCR = pd.read_csv(f"{TOP}/a9_snrna_screen.tsv", sep="\t").set_index("gene")
CONC = pd.read_csv(f"{TOP}/a9_concordance_screened.tsv", sep="\t")
OX = pd.read_csv(f"{TOP}/a2_oxphos_level.tsv", sep="\t")
RF = pd.read_csv(f"{TOP}/a0_primary_refit_summary.tsv", sep="\t").set_index(["refit", "gene_set"])
PROG = [l.strip() for l in open(f"{OUT}/geneset_primary_program_qcpass_down.txt") if l.strip() and not l.startswith("#")]
N_PROG = len(PROG)
assert set(PROG) == set(H[H.qc_pass & (H.logFC < 0)].index), "score gene set != QC-passing lower-in-CS hits"
CH = {}
for tag in ("DCM", "HCM"):
    t = pd.read_csv(f"data/published/chaffin2022_cardiomyocyte_{tag}vsNF.tsv", sep="\t")
    CH[tag] = t[t.bg_contamination_flag == 0].set_index("gene")


def qcv(univ, d, stat, col="qc_passing_only"):
    r = qc[(qc.analysis == "cross_platform_snrna") & (qc.gene_universe == univ) & (qc.direction == d) &
           (qc.statistic == stat)]
    assert len(r) == 1
    return float(r[col].iloc[0])


def concv(screen, d, gs, col):
    r = CONC[(CONC.screen == screen) & (CONC.direction == d) & (CONC.gene_set == gs)]
    assert len(r) == 1
    return float(r[col].iloc[0])


# ------------------------------------------------------------------ panel helpers
def panel_A(ax, sec, img_path, scale, title, spot_diam):
    im = Image.open(img_path).convert("RGB")
    s = sp[sp.section == sec]
    x, y = s.px_x.values * scale, s.px_y.values * scale
    pad = 40
    x0, x1 = max(0, x.min() - pad), min(im.size[0], x.max() + pad)
    y0, y1 = max(0, y.min() - pad), min(im.size[1], y.max() + pad)
    ax.imshow(im, extent=(0, im.size[0], im.size[1], 0))
    cat = np.where(s.immune_dominant, "imm", np.where(s.cm_dominant, "cm", "oth"))
    d = spot_diam * scale
    for k, c in (("cm", C_SEL), ("imm", C_IMM)):  # other-lineage spots left unmarked
        m = cat == k
        ec = EllipseCollection(np.full(m.sum(), d), np.full(m.sum(), d), np.zeros(m.sum()), units="xy",
                               offsets=np.c_[x[m], y[m]], offset_transform=ax.transData, facecolors=c,
                               edgecolors="none", alpha=0.85 if k != "oth" else 0.55)
        ec.set_rasterized(True); ax.add_collection(ec)
    ax.set_xlim(x0, x1); ax.set_ylim(y1, y0); ax.set_xticks([]); ax.set_yticks([])
    for sp_ in ax.spines.values(): sp_.set_visible(False)
    # 1 mm scale bar
    ppmm = float(s.px_per_um.iloc[0]) * 1000 * scale
    bx, by = x1 - ppmm - 0.04 * (x1 - x0), y1 - 0.04 * (y1 - y0)
    ax.plot([bx, bx + ppmm], [by, by], color="k", lw=1.2)
    ax.text(bx + ppmm / 2, by - 0.015 * (y1 - y0), "1 mm", ha="center", va="bottom", fontsize=6)
    ax.set_title(title, loc="left", pad=2)
    n = dict(pd.Series(cat).value_counts())
    return {k: int(n.get(k, 0)) for k in ("cm", "imm", "oth")}


SC_LAB = {"score_program": f"Score of genes\nlower in CS (n = {N_PROG})", "score_ACADVL": it("ACADVL")}
SLOPE_KEY = {"score_program": "program", "score_ACADVL": "ACADVL"}   # row names in a2_perpatient_slopes.tsv
ORDER = list(bins.bin)


def panel_B(ax, score, show_x, slope_text=True, fs=5.3, headroom=0.06):
    xs = np.arange(len(ORDER))
    ref = cref.loc[score]
    ax.axhspan(ref["min"], ref["max"], color=C_CMP, alpha=0.15, lw=0)
    ax.axhline(ref["mean_across_patients"], color=C_CMP, lw=0.9)
    for p, d in bpp.groupby("patient"):
        xi = [ORDER.index(b) for b in d.bin]
        ax.plot(xi, d[score], color=C_CS, lw=0.5, alpha=0.35, marker="o", ms=1.5, mew=0)
    ok = bins[f"{score}_mean"].notna().values
    m, lo, hi = (bins[f"{score}_{k}"].values for k in ("mean", "ci_low", "ci_high"))
    assert ok.all(), "unpopulated bin reached the plot"
    ax.fill_between(xs, lo, hi, color=C_CS, alpha=0.18, lw=0)
    ax.plot(xs, m, color=C_CS, lw=1.6, marker="s", ms=3)
    lo_all = np.nanmin([np.nanmin(bpp[score]), ref["min"]]); hi_all = np.nanmax([np.nanmax(bpp[score]), ref["max"]])
    ax.set_ylim(lo_all - 0.06 * (hi_all - lo_all), hi_all + headroom * (hi_all - lo_all))
    ax.set_xticks(xs); ax.set_xticklabels(ORDER if show_x else [])
    ax.set_xlim(-0.5, len(ORDER) - 0.5)
    ax.set_ylabel(SC_LAB[score] + "\n(log1p CP10K)")
    r = slp.loc[SLOPE_KEY[score]]
    n, kneg = int(r.n_patients), int(r.n_negative)
    if slope_text:   # above the plotting area, off the data
        ax.text(1.0, 1.01, f"Per-patient slope {r.mean_slope_per_mm:+.3f} per mm "
                           f"(95% CI {r.ci_low:+.3f} to {r.ci_high:+.3f}; {kneg}/{n} negative)",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=fs)
    for i, b in enumerate(ORDER):
        val("B", f"{score}|{b}|pooled_mean", m[i])
    val("B", f"{score}|comparator_mean", ref["mean_across_patients"])
    val("B", f"{score}|perpatient_slope", r.mean_slope_per_mm); val("B", f"{score}|perpatient_slope_ci_low", r.ci_low)
    val("B", f"{score}|perpatient_slope_ci_high", r.ci_high); val("B", f"{score}|perpatient_n", n)
    val("B", f"{score}|perpatient_n_negative", kneg)


def counts_row(ax, fs=5.5):
    tr = mpl.transforms.blended_transform_factory(ax.transData, ax.transAxes)
    for i, r in enumerate(bins.itertuples()):
        ax.text(i, -0.17, f"{int(r.n_spots)}/{int(r.n_patients)}", transform=tr, ha="center", va="top", fontsize=fs, color="0.3")
    ax.text(-0.55, -0.17, "spots/patients", transform=tr, ha="right", va="top", fontsize=fs, color="0.3")


# ------------------------------------------------------------------ Panel C / D inputs (screen 3a applied)
q = H[H.qc_pass].copy()
q = q[q.index.isin(SN.index)]
q["sn"] = SN.loc[q.index, "lfc_pool"]; q["fdr"] = SN.loc[q.index, "fdr_pool"]
scr_q = SCR[SCR.qc_pass & SCR.testable_snrna]
assert sorted(scr_q.index) == sorted(q.index), "screen table and plotted QC-passing set differ"
assert np.allclose(scr_q.loc[q.index, "lfc_pool"], q.sn)
q["conc"] = np.sign(q.sn) == np.sign(q.logFC)
q["rep"] = q.conc & (q.fdr < 0.05)
q["fail"] = q.sn.abs() > SCREEN
assert (q.fail == SCR.loc[q.index, "fail_3a"].astype(bool)).all()
q["oor"] = q.sn.abs() > YLIM_C
assert (q.oor == q.fail).all(), "a gene beyond the axis passes the screen, or a screened gene is in range"
q["y"] = q.sn.clip(-YLIM_C, YLIM_C)
q["klass"] = np.where(q.rep, "rep", np.where(q.conc, "same", "opp"))
qs = q[~q.fail]
AG = {}
for d, m in (("down", qs.logFC < 0), ("up", qs.logFC > 0)):
    sub = qs[m]
    n_same = int(sub.conc.sum())
    assert len(sub) == concv("screened_3a", d, "qc_passing_only", "n_testable")
    assert n_same == concv("screened_3a", d, "qc_passing_only", "n_same_sign")
    full = q[(q.logFC < 0) if d == "down" else (q.logFC > 0)]   # method check: unscreened rho reproduces A3
    assert np.isclose(stats.spearmanr(full.logFC, full.sn).statistic, qcv("all_fdr5_hits", d, "spearman_rho"), atol=1e-9)
    AG[d] = dict(frac=n_same / len(sub), bg=concv("screened_3a", d, "qc_passing_only", "background_frac"),
                 rho=stats.spearmanr(sub.logFC, sub.sn).statistic, n=len(sub),
                 frac_un=concv("unscreened", d, "qc_passing_only", "frac_same_sign"),
                 bg_un=concv("unscreened", d, "qc_passing_only", "background_frac"))
CLS = (("opp", C_BG, 7, 0, "Opposite sign"), ("same", "#6BAED6", 7, 0, "Same sign"),
       ("rep", C_CMP, 12, 0.3, "Same sign, snRNA FDR <5%"))

genesD = list(SCR[SCR.in_concordant_core_v2_screened].sort_values("spatial_logFC").index)  # user decision: core after 3a
_core = pd.read_csv(f"{TOP}/primary_v2/q1_core_snrna_pergene.tsv", sep="\t").set_index("gene")
_core18 = set(_core.index[_core.replicates_fdr5])
assert set(genesD) == {g for g in _core18 if not bool(SCR.loc[g, "fail_3a"])} and set(genesD) <= set(H.index)
val("D", "n_core_removed_screen_3a", len(_core18) - len(genesD))
# snRNA vs NF column: no plotted value beyond the 3a magnitude (user check, v5)
_nfmax = max(abs(NF.loc[g, "logFC"]) for g in genesD if g in NF.index)
assert _nfmax <= SCREEN, f"snRNA vs NF |log2FC| {_nfmax:.2f} beyond {SCREEN}"
val("D", "max_abs_snrna_vs_nf_logFC", _nfmax)
val("D", "min_snrna_vs_nf_AveExpr", min(NF.loc[g, "AveExpr"] for g in genesD if g in NF.index))


def scatter_C(ax, qq, labels, extra_label=None, fs=5.3, s_mult=1.0):
    ax.axhline(0, color="0.7", lw=0.5); ax.axvline(0, color="0.7", lw=0.5)
    for k, c, sz, lw, lb in CLS:
        m = qq.klass == k
        mi = m & ~qq.oor
        ax.scatter(qq.logFC[mi], qq.y[mi], s=sz * s_mult, c=c, lw=lw, ec="k" if lw else "none", label=lb)
    for i, (g, r) in enumerate(qq[qq.oor].iterrows()):   # screened (3a): open neutral outward triangle, no category
        ax.scatter([r.logFC], [r.y], marker="^" if r.sn > 0 else "v", s=16 * s_mult, facecolors="none",
                   edgecolors=C_SCR, lw=0.6, clip_on=False, zorder=5,
                   label=f"snRNA |log2FC| >{SCREEN:g}, screened" if i == 0 else None)
    ax.set_ylim(-YLIM_C, YLIM_C)
    dx = 0.02 * np.ptp(qq.logFC.values)
    for sgn in (1, -1):   # off-scale: value only (no gene names); two-row stagger, alternating sides within clusters
        o = qq[qq.oor & (np.sign(qq.sn) == sgn)].sort_values("logFC")
        cl = np.cumsum(np.r_[0, np.diff(o.logFC.values) > 5 * dx]) if len(o) else []
        k_in = 0
        for i, (g, r) in enumerate(o.iterrows()):
            k_in = 0 if i == 0 or cl[i] != cl[i - 1] else k_in + 1
            ax.text(r.logFC + (-dx if k_in % 2 == 0 else dx), sgn * (YLIM_C - 0.2 - 0.22 * (cl[i] % 2)), f"{r.sn:+.1f}",
                    fontsize=fs, ha="right" if k_in % 2 == 0 else "left", va="center")
    txt = [ax.text(qq.loc[g, "logFC"], qq.loc[g, "y"], g, fontsize=fs, style="italic") for g in labels]
    if extra_label is not None:
        txt.append(ax.text(qq.loc[extra_label, "logFC"], qq.loc[extra_label, "y"], extra_label, fontsize=fs + 0.7, style="italic"))
    adjust_text(txt, x=qq.logFC.values, y=qq.y.values, ax=ax, arrowprops=dict(arrowstyle="-", color="0.4", lw=0.4))


# ------------------------------------------------------------------ figure (2 x 2)
fig = plt.figure(figsize=(180 / 25.4, 190 / 25.4))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.22], hspace=0.3, wspace=0.34,
                      left=0.085, right=0.985, top=0.955, bottom=0.075)
gA = gs[0, 0].subgridspec(1, 2, wspace=0.04)
axA1, axA2 = fig.add_subplot(gA[0]), fig.add_subplot(gA[1])
cs_scale = imeta[CS_SEC]["hires_scalef"] * imeta[CS_SEC]["extra_downscale"]
nA_cs = panel_A(axA1, CS_SEC, f"data/eschf/img/{CS_SEC}.jpg", cs_scale, "Cardiac sarcoidosis",
                 imeta[CS_SEC]["spot_diameter_fullres"])
cmp_scale = imeta[COMP_SEC]["hires_scalef"] * imeta[COMP_SEC]["extra_downscale"]
nA_cmp = panel_A(axA2, COMP_SEC, f"data/eschf/img/{COMP_SEC}.jpg", cmp_scale, "Comparator (ARVC)",
                 imeta[COMP_SEC]["spot_diameter_fullres"])
for k, v in nA_cs.items(): val("A", f"{CS_SEC}|{k}", v)
for k, v in nA_cmp.items(): val("A", f"{COMP_SEC}|{k}", v)
letter(axA1, "A", dx=-0.04)
axA1.legend(handles=[Line2D([], [], ls="", marker="o", ms=4, mfc=c, mec="none", label=l) for c, l in (
    (C_SEL, "Cardiomyocyte-dominant (analysed)"), (C_IMM, "Immune-dominant"))] +
    [Line2D([], [], ls="", label="Unmarked: other dominant lineage")],
    loc="upper left", bbox_to_anchor=(0.0, -0.01), ncol=2, frameon=False, handletextpad=0.2, columnspacing=0.8,
    fontsize=5.3)

gB = gs[0, 1].subgridspec(2, 1, hspace=0.3)
axB1, axB2 = fig.add_subplot(gB[0]), fig.add_subplot(gB[1])
panel_B(axB1, "score_program", False, headroom=0.3)
panel_B(axB2, "score_ACADVL", True)
axB1.set_title("Expression by distance from the nearest immune-dominant spot", loc="left", pad=11)
counts_row(axB2)
axB2.set_xlabel("Distance to nearest immune-dominant spot (mm)", labelpad=9)
axB1.legend(handles=[Line2D([], [], color=C_CS, lw=1.6, marker="s", ms=3, label="CS mean of patients (95% CI)"),
                     Line2D([], [], color=C_CS, lw=0.5, alpha=0.5, label="CS patient"),
                     Patch(fc=C_CMP, alpha=0.3, ec=C_CMP, label="Comparator mean (range)")],
            loc="upper left", frameon=False, fontsize=5.3, ncol=3, bbox_to_anchor=(0.0, 1.0), handlelength=1.6,
            columnspacing=1.0, handletextpad=0.4)   # single row in headroom above the comparator band
letter(axB1, "B", dx=-0.15, dy=1.1)

axC = fig.add_subplot(gs[1, 0])
labC = [g for g in genesD if g in qs.index]
scatter_C(axC, q, labC)
axC.set_xlabel("Spatial log2 fold change, CS vs comparators")
axC.set_ylabel("snRNA pooled log2 fold change\n(CS vs DCM and ARVC)")
axC.text(0.01, 0.02, f"Lower in CS: {AG['down']['frac']:.0%} same sign (background {AG['down']['bg']:.0%}), ρ = {AG['down']['rho']:+.2f}\n"
                     f"Higher in CS: {AG['up']['frac']:.0%} same sign (background {AG['up']['bg']:.0%}), ρ = {AG['up']['rho']:+.2f}",
         transform=axC.transAxes, ha="left", va="bottom", fontsize=5.3)
axC.legend(loc="upper left", bbox_to_anchor=(0.0, 0.86), frameon=False, fontsize=5.3, handletextpad=0.2)
axC.set_title("Spatial vs snRNA effects, QC-passing hits", loc="left", pad=8)
axC.margins(x=0.06)
letter(axC, "C", dx=-0.13, dy=1.05)
val("C", "n_qcpass_genes_plotted", len(q)); val("C", "n_labelled", len(labC)); val("C", "n_replicating", int(q.rep.sum()))
val("C", "n_out_of_range", int(q.oor.sum())); val("C", "n_screen_3a_fail", int(q.fail.sum()))
for d in ("down", "up"):
    val("C", f"screened|frac_same_sign_{d}", AG[d]["frac"]); val("C", f"screened|background_{d}", AG[d]["bg"])
    val("C", f"screened|rho_{d}", AG[d]["rho"]); val("C", f"screened|n_{d}", AG[d]["n"])
    val("C", f"unscreened|frac_same_sign_{d}", AG[d]["frac_un"]); val("C", f"unscreened|background_{d}", AG[d]["bg_un"])
for g, r in q[q.oor].iterrows():
    val("C", f"out_of_range|{g}|snrna_lfc", r.sn)
    val("C", f"out_of_range|{g}|drawn_open_neutral", int(r.fail))

# D: dot matrix, concordant genes of the section-stable core after screen 3a
gD = gs[1, 1].subgridspec(1, 2, width_ratios=[1, 0.5], wspace=0.02)
axD = fig.add_subplot(gD[0]); axDL = fig.add_subplot(gD[1]); axDL.axis("off")
cols = [("Spatial\nCS vs\ncomp.", lambda g: (H.loc[g, "logFC"], H.loc[g, "adj.P.Val"])),
        ("snRNA\nvs DCM\n+ARVC", lambda g: (SN.loc[g, "lfc_pool"], SN.loc[g, "fdr_pool"])),
        ("snRNA\nvs NF", lambda g: (NF.loc[g, "logFC"], NF.loc[g, "adj.P.Val"]) if g in NF.index else None),
        ("Chaffin\nDCM\nvs NF", lambda g: (CH["DCM"].loc[g, "logFC"], CH["DCM"].loc[g, "adj_p"]) if g in CH["DCM"].index else None),
        ("Chaffin\nHCM\nvs NF", lambda g: (CH["HCM"].loc[g, "logFC"], CH["HCM"].loc[g, "adj_p"]) if g in CH["HCM"].index else None)]
for j, (cname, f) in enumerate(cols):
    for i, g in enumerate(genesD):
        r = f(g)
        if r is None:
            axD.scatter(j, i, marker="s", s=22, facecolors="none", edgecolors="0.5", lw=0.6)
            val("D", f"{g}|{cname.replace(chr(10), ' ')}|logFC", np.nan)
            continue
        lfc, p = r
        axD.scatter(j, i, s=6 + 34 * min(-np.log10(max(p, 1e-10)), 10) / 10, c=C_LO if lfc < 0 else C_HI, lw=0)
        val("D", f"{g}|{cname.replace(chr(10), ' ')}|logFC", lfc)
SYM = {}
for g in genesD:
    r = SCR.loc[g]
    dag, ddag, sec = (r.flag_locus_quant_outlier is True or r.flag_locus_quant_outlier == True), (not bool(r.in_qc_universe)), (r.flag_segdup == True)
    assert not (dag and ddag) and not (sec and ddag)
    assert not bool(SCR.loc[g, "fail_3a"])
    SYM[g] = ("†" if dag else "") + ("‡" if ddag else "") + ("§" if sec else "")
    val("D", f"{g}|flag_dagger_locus_outlier", int(dag)); val("D", f"{g}|flag_ddagger_out_of_universe", int(ddag))
    val("D", f"{g}|flag_section_segdup", int(sec))
used = set("".join(SYM.values()))
axD.set_yticks(range(len(genesD))); axD.set_yticklabels([f"{g} {SYM[g]}".strip() for g in genesD], style="italic", fontsize=5.5)
axD.set_xticks(range(len(cols))); axD.set_xticklabels([c for c, _ in cols], fontsize=5.2)
axD.set_xlim(-0.6, len(cols) - 0.4); axD.set_ylim(len(genesD) - 0.5, -0.5)
axD.spines["left"].set_visible(False); axD.spines["bottom"].set_visible(False); axD.tick_params(length=0)
axD.set_title(f"{len(genesD)} section-stable concordant genes", loc="left", pad=2)
val("D", "n_concordant_genes", len(genesD))
for s_ in "†‡§": val("D", f"n_symbol|{s_}", sum(s_ in v for v in SYM.values()))
hd = [Line2D([], [], ls="", marker="o", ms=4, mfc=C_LO, mec="none", label="Lower in disease"),
      Line2D([], [], ls="", marker="o", ms=4, mfc=C_HI, mec="none", label="Higher in disease"),
      Line2D([], [], ls="", marker="s", ms=4, mfc="none", mec="0.5", label="Not tested, or not in\npublished DEG list"),
      Line2D([], [], ls="", marker="o", ms=2, mfc="0.3", mec="none", label="FDR 0.1"),
      Line2D([], [], ls="", marker="o", ms=5.5, mfc="0.3", mec="none", label="FDR ≤1e-10")]
if "†" in used: hd.append(Line2D([], [], ls="", label="† locus-quantification\n   outlier flag"))
if "‡" in used: hd.append(Line2D([], [], ls="", label="‡ outside locus-QC\n   universe (not evaluated)"))
if "§" in used: hd.append(Line2D([], [], ls="", label="§ segmental-duplication\n   flag"))
axDL.legend(handles=hd, loc="upper left", bbox_to_anchor=(0.0, 1.0), frameon=False, fontsize=5.2,
            handletextpad=0.3, labelspacing=0.7, handlelength=1.0)
letter(axD, "D", dx=-0.3, dy=1.02)

fig.savefig(f"{FIG}/eschf_figure1.pdf", dpi=600)
fig.savefig(f"{FIG}/eschf_figure1.png", dpi=300)

# ------------------------------------------------------------------ graphical abstract (11 x 18 cm)
ga = plt.figure(figsize=(11 / 2.54, 18 / 2.54))
g2 = ga.add_gridspec(2, 1, height_ratios=[2, 1.25], hspace=0.3, left=0.2, right=0.96, top=0.92, bottom=0.07)
g2t = g2[0].subgridspec(2, 1, hspace=0.12)
ga.suptitle("Cardiomyocyte expression in cardiac sarcoidosis", fontsize=8, x=0.06, ha="left", y=0.985)
a1, a2, a3 = ga.add_subplot(g2t[0]), ga.add_subplot(g2t[1]), ga.add_subplot(g2[1])
panel_B(a1, "score_program", False, slope_text=False, headroom=0.3); panel_B(a2, "score_ACADVL", True, slope_text=False)
a1.set_title("Cardiomyocyte scores by distance from immune-dominant spots", loc="left", pad=2)
counts_row(a2, fs=5.0)
a2.set_xlabel("Distance to nearest immune-dominant spot (mm)", labelpad=13)
a1.legend(handles=[Line2D([], [], color=C_CS, lw=1.6, marker="s", ms=3, label="Sarcoidosis"),
                   Patch(fc=C_CMP, alpha=0.3, ec=C_CMP, label="Other cardiomyopathies")],
          loc="upper left", bbox_to_anchor=(0.0, 1.0), frameon=False, fontsize=6, ncol=2, columnspacing=1.2)
ga_lab = "ACADVL" if "ACADVL" in qs.index else None   # the only gene label allowed
scatter_C(a3, q, [], extra_label=ga_lab, fs=5.5, s_mult=0.9)
a3.legend(loc="upper left", bbox_to_anchor=(0.0, 0.86), frameon=False, fontsize=6, handletextpad=0.2)
a3.set_xlabel("Spatial log2 fold change"); a3.set_ylabel("snRNA log2 fold change")
a3.set_title("Spatial vs snRNA effects, QC-passing hits", loc="left", pad=8)
a3.text(0.01, 0.02, f"Lower-in-CS genes: {AG['down']['frac']:.0%} same sign (background {AG['down']['bg']:.0%})",
        transform=a3.transAxes, ha="left", va="bottom", fontsize=6)
ga.savefig(f"{FIG}/eschf_graphical_abstract.png", dpi=600)

# ------------------------------------------------------------------ legend + alt text (numbers from the data above)
mx = float(pd.read_csv(f"{OUT}/a2_model_input.tsv.gz", sep="\t").dist_gran_mm.max())
val("B", "max_distance_cm_spot_mm", mx); val("B", "n_score_genes", N_PROG)
nB_pat = int(slp.loc["program", "n_patients"])
pr = pmap[pmap.retained == "yes"]   # sections in the primary fit (CS_109: no CM-dominant spots)
ncs_fit = int((pr.arm == "CS").sum()); ncs_pat_fit = int(pr[pr.arm == "CS"].patient.nunique())
ncmp_fit = int((pr.arm == "COMPARATOR").sum()); ncmp_pat_fit = int(pr[pr.arm == "COMPARATOR"].patient.nunique())
n_cs_excl = int(((pmap.arm == "CS") & (pmap.retained != "yes")).sum())
val("legend", "n_cs_sections_fit", ncs_fit); val("legend", "n_cs_patients_fit", ncs_pat_fit)
val("legend", "n_comparator_sections_fit", ncmp_fit); val("legend", "n_comparator_patients_fit", ncmp_pat_fit)
pp = RF.loc[("patient_pseudobulk", "all_tested")]
val("legend", "patient_pseudobulk|n_section_hits", pp.n_section_hits)
val("legend", "patient_pseudobulk|section_hits_still_fdr5", pp.section_hits_still_fdr5)
val("legend", "patient_pseudobulk|sign_concordant_section_hits", pp.sign_concordant_section_hits)
# B legend pathway sentence (user wording, v4): KEGG TCA / FA degradation / OXPHOS not lower in CS once fit hits are
# excluded. Asserted against A14: for minus_primary and minus_union, no 95% CI lies entirely below zero.
A14 = pd.read_csv(f"{TOP}/a14_pathway_hits_removed.tsv", sep="\t")
_pw = {"tca_cycle_hsa00020", "fatty_acid_degradation_hsa00071", "oxphos_hsa00190"}
_hr = A14[A14.version.isin(["minus_primary", "minus_union"])]
assert set(_hr.set) == _pw and len(_hr) == 6 and (~_hr.too_small_to_interpret.astype(bool)).all()
assert (_hr.diff_ci_high > 0).all(), "pathway sentence ('not lower in CS once fit hits were excluded') no longer true"
for _, r in _hr.iterrows():
    for k in ("diff_cs_minus_comparator", "diff_ci_low", "diff_ci_high"):
        val("legend", f"pathway_hits_removed|{r.set}|{r.version}|{k}", r[k])
n_loo = len(glob.glob(f"{TOP}/primary_v2/loo/q1_raw_drop_*.tsv"))
assert n_loo == 6
val("legend", "n_section_refits_core", n_loo)
fails = q[q.fail].sort_values("sn")
k_fail = len(fails)
# C legend source-absence sentence (user wording, v5): every 3a failure carries the 3b flag (one study near-absent,
# others >= 2 log2 CPM higher), as for the FPGT reference gene of the FPGT-TNNI3K locus.
_f3 = SCR[SCR.fail_3a.astype(bool)]
assert len(_f3) == k_fail and _f3.flag_3b.astype(bool).all(), "not all screened genes carry the 3b source-absence flag"
_SA = pd.read_csv(f"{TOP}/a9_source_absence.tsv", sep="\t").set_index("gene")
assert bool(_SA.loc["FPGT", "flag_3b"]), "FPGT reference not flagged by 3b"
for g, r in _f3.iterrows():
    val("legend", f"source_absence|{g}|flag_3b", int(bool(r.flag_3b)))
    val("legend", f"source_absence|{g}|absent_mean_log2cpm", r.absent_mean)
    val("legend", f"source_absence|{g}|others_mean_log2cpm", r.others_mean)
val("legend", "source_absence|FPGT|flag_3b", int(bool(_SA.loc["FPGT", "flag_3b"])))
_NUM = {0: "No", 1: "One", 2: "Two", 3: "Three", 4: "Four", 5: "Five", 6: "Six", 7: "Seven", 8: "Eight", 9: "Nine"}
n_fail_word = ("All " + _NUM[k_fail].lower()) if k_fail > 1 else _NUM[k_fail]
fail_txt = ", ".join(f"*{g}* ({r.sn:+.1f})" for g, r in fails.iterrows())
dag_g = [g for g in genesD if "†" in SYM[g]]; ddag_g = [g for g in genesD if "‡" in SYM[g]]; sec_g = [g for g in genesD if "§" in SYM[g]]
fmt = lambda L: ", ".join(f"*{g}*" for g in L) if L else "none"
pB = slp.loc["program"]; aB = slp.loc["ACADVL"]
A, U = AG["down"], AG["up"]
legend = f"""# Figure 1. Cardiomyocyte expression in cardiac sarcoidosis by distance from immune-dominant spots and across platforms

All panels use the primary contrast: cardiac sarcoidosis (CS) against comparator cardiomyopathies over **all**
cardiomyocyte (CM)-dominant Visium spots, with no exclusion by distance to immune-dominant spots ({ncs_fit} CS
sections from {ncs_pat_fit} patients; {ncmp_fit} comparator sections from {ncmp_pat_fit} patients; section-level
limma-voom; {n_cs_excl} CS section with no CM-dominant spot passing the UMI floor is not included). At patient level
(sections summed per patient; {ncs_pat_fit} CS and {ncmp_pat_fit} comparator patients), {int(pp.section_hits_still_fdr5)} of the
{int(pp.n_section_hits)} section-level hits remained at FDR <5% and {int(pp.sign_concordant_section_hits)} of {int(pp.n_section_hits)} kept their sign. The earlier fit restricted
to CM-dominant spots ≥0.5 mm from an immune-dominant spot is reported as a sensitivity analysis.

**(A)** H&E of one CS section ({CS_SEC}) and one comparator section ({COMP_SEC}, arrhythmogenic cardiomyopathy),
with Visium spots that pass the 2000-UMI floor coloured by dominant lineage. Green, CM-dominant spots (all analysed;
{nA_cs['cm']} in {CS_SEC}, {nA_cmp['cm']} in {COMP_SEC}); orange, immune-dominant spots ({nA_cs['imm']} in
{CS_SEC}; {nA_cmp['imm']} in {COMP_SEC}). Spots of other dominant lineages are unmarked. Scale bars, 1 mm.
**(B)** Score of genes lower in CS (top; n = {N_PROG}) and *ACADVL* (bottom) in CM-dominant CS spots, by distance to
the nearest immune-dominant spot. The score averages the {N_PROG} QC-passing genes lower in CS in the primary fit; its
level difference from comparators reflects how the genes were selected, and the panel shows how it varies with
distance. Score, mean log1p CP10K over member genes on depth-matched counts. Thin lines, per-patient bin means; thick
line and band, mean across patients and 95% t-interval. Blue line and band, comparator mean and range of per-patient
means. Counts beneath the axis are spots/patients per bin. Only populated bins are shown. The maximum distance from
any CS CM-dominant spot to the nearest immune-dominant spot was {mx:.2f} mm; no CS section contained CM-dominant
myocardium beyond 2 mm from an immune-dominant spot. Text above each plot, per-patient slope: for each CS patient
with CM-dominant spots in at least two populated bins (n = {nB_pat}), ordinary least squares on that patient's spots,
outcome ~ distance + CM-structural score + z(log10 UMI), with a section indicator for patients with more than one
section; shown are the mean of the per-patient slopes, its 95% t-interval and the number of negative slopes (score
{pB.mean_slope_per_mm:+.3f} per mm, 95% CI {pB.ci_low:+.3f} to {pB.ci_high:+.3f}; *ACADVL* {aB.mean_slope_per_mm:+.3f} per mm,
95% CI {aB.ci_low:+.3f} to {aB.ci_high:+.3f}). Pathway scores defined independently of the fit (KEGG TCA cycle, fatty
acid degradation and oxidative phosphorylation) were not lower in CS once fit hits were excluded.
**(C)** Spatial log2 fold change (CS vs comparators) against pooled snRNA log2 fold change (CS vs DCM and CS vs ARVC,
inverse-variance pooled) for spatial hits (FDR <5%) passing locus QC and testable in snRNA (n = {len(q)}). Because the
CS and comparator single-nucleus cohorts come from different studies, pooled single-nucleus effects beyond ±{SCREEN:g}
log2 were treated as potentially affected by cross-study quantification. Genes with pooled single-nucleus effects beyond ±{SCREEN:g} log2 (n = {k_fail}) are shown at the axis
edge and excluded from the agreement fractions: {fail_txt}. {n_fail_word} screened genes were near-absent from one
source study's quantification while present in the others, as previously observed at the FPGT-TNNI3K locus; the
±{SCREEN:g} log2 threshold is a pre-stated sensitivity screen, not a biological limit. Triangles point outward, are drawn open in grey without a sign
category and are labelled with their value only; no other gene lies beyond the ±{YLIM_C:g} y-axis. Labelled, the genes of
(D) that pass locus QC and the screen. Text, fraction of genes with the same sign in each direction after the screen, beside the background fraction
from the label-permutation null (lower in CS {A['frac']:.0%}, n = {A['n']}, vs {A['bg']:.0%}; higher in CS {U['frac']:.0%}, n = {U['n']},
vs {U['bg']:.0%}; without the screen {A['frac_un']:.0%} vs {A['bg_un']:.0%} and {U['frac_un']:.0%} vs {U['bg_un']:.0%}); direction agreement is reported
descriptively. ρ, Spearman correlation of fold changes after the screen.
**(D)** Concordant genes of the section-stable core; genes failing locus QC are marked. Shown
across contrasts are the {len(genesD)} genes of the section-stable core (spatial FDR <5% in the full fit and in each of the
{n_loo} refits leaving out one comparator section) with the same sign and snRNA FDR <5% in the pooled single-nucleus
contrast, after the screen in (C) ({fmt(sorted(_core18 - set(genesD)))} removed). Colour, direction of log2 fold change; dot size, −log10 FDR (capped at 10). Open squares, gene not
tested, or absent from the published differential-gene list (Chaffin et al. 2022, cardiomyocytes; genes with the
background-contamination flag excluded). †, locus-quantification outlier flag ({fmt(dag_g)}); ‡, gene outside the
locus-QC universe, not evaluated by locus QC ({fmt(ddag_g)}).{'' if not sec_g else ' §, segmental-duplication flag (' + fmt(sec_g) + ').'}

## Alt text
Four-panel figure. A: two H&E sections with coloured spots marking analysed cardiomyocyte-dominant spots and
immune-dominant spots. B: two line plots, a score of genes lower in sarcoidosis and ACADVL, against distance from the
nearest immune-dominant spot over {len(ORDER)} distance bins up to 2 mm, with sarcoidosis patient means and the comparator
mean and range. C: scatter of spatial against single-nucleus fold changes on a limited vertical axis, with {k_fail}
screened genes marked by open grey triangles at the edge. D: dot matrix of {len(genesD)} section-stable concordant genes across five
contrasts, with many open squares in the published-reference columns.

## Graphical abstract (11 × 18 cm, width × height)
Panels B (score of genes lower in CS, *ACADVL*) and C (sign agreement after the snRNA screen) of Figure 1, redrawn with
the same populated bins and off-scale rule; off-scale points carry their value only, and *ACADVL* is the only gene
labelled.
"""
open(f"{FIG}/eschf_figure1_legend.md", "w").write(legend)
pd.DataFrame(VALS).to_csv(f"{TOP}/a7_panel_values.tsv", sep="\t", index=False)
print(len(VALS), "values;", "A:", nA_cs, nA_cmp, "| D symbols:", {g: s for g, s in SYM.items() if s})
print("C triangles:", q[q.oor][["logFC", "sn", "klass"]].round(2).to_dict("index"))
print("C labels:", len(labC), "| agreement:", {d: {k: round(v, 3) for k, v in AG[d].items()} for d in AG})
