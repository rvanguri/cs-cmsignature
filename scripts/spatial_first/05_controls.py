#!/usr/bin/env python3
"""Technical controls on the spatial Q1 contrast.

Each control asks whether the CS-vs-comparator difference could be an artefact of the
way the two cohorts were generated rather than a difference between the tissues.

  panelB          CS vs only the comparator sections run on probe panel B
  panelmatched    CS vs only the comparator sections run on the same panel as CS
  purityadj       CS vs all comparators, adjusting for per-section cardiomyocyte purity
  panel_contrast  panel B vs panel A among comparators only, i.e. the size of the
                  panel effect measured directly, with disease held out of the contrast
  noMT            mitochondrially encoded genes removed and TMM normalisation recomputed
                  from scratch, so a shift in mitochondrial content cannot drive the fit

All four reuse the committed pseudobulk matrix and the detectability table; the first
four call 01b_limma_twogroup.R, noMT calls 01_limma_q1q2.R on an MT-stripped matrix.

Outputs (results/spatial_first/)
  q1_panelB.tsv, q1_panelmatched.tsv, q1_purityadj.tsv, ctrl_panelB_vs_panelA.tsv
  q1_raw_T2000_d25_noMT.tsv, q1_raw_T2000_d25_noMTrenorm.tsv
  q1_noMT_summary.tsv       tested / significant / direction with and without MT genes
  q1_noMT_comparison.tsv    per-gene logFC and FDR, original vs MT-free
  q1_control_comparison.tsv per-gene logFC and FDR across the three Q1 control arms,
                            restricted to genes significant in the main fit
"""
import argparse
import os
import subprocess
import tempfile

import numpy as np
import pandas as pd

DATA = "data/spatial_first"
RES = "results/spatial_first"
SCR = "scripts/spatial_first"
PB = f"{DATA}/pseudobulk_T2000.tsv.gz"
DET = f"{DATA}/detection_fraction_T2000.tsv.gz"
THRESH = "0.25"

# name -> (meta file, group 1, group 2, detectability columns, covariate)
TWOGROUP = {
    "q1_panelB":             ("meta_q1_panelB.tsv", "CS", "COMPB", "CS,COMPARATOR", None),
    "q1_panelmatched":       ("meta_q1_panelmatched.tsv", "CS", "COMPA", "CS,COMPARATOR", None),
    "q1_purityadj":          ("meta_q1_purity.tsv", "CS", "COMPARATOR", "CS,COMPARATOR", "cm_purity"),
    # comparators only, so detectability is assessed in the comparator arm alone
    "ctrl_panelB_vs_panelA": ("meta_panel.tsv", "panelB", "panelA", "COMPARATOR", None),
}


# the R stages usually live in a separate conda environment
RSCRIPT = os.environ.get("RSCRIPT", "Rscript")


def run(cmd):
    subprocess.run(cmd, shell=True, check=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-fits", action="store_true",
                    help="only rebuild the comparison/summary tables from existing fits")
    args = ap.parse_args()

    if not args.skip_fits:
        for name, (meta, g1, g2, detcols, covar) in TWOGROUP.items():
            cov = f" --covar {covar}" if covar else ""
            run(f"{RSCRIPT} {SCR}/01b_limma_twogroup.R --pb {PB} --det {DET} "
                f"--meta {DATA}/{meta} --thresh {THRESH} --g1 {g1} --g2 {g2} "
                f"--detcols {detcols}{cov} --outfile {RES}/{name}.tsv")

        # --- MT-free refit --------------------------------------------------------
        pb = pd.read_csv(PB, sep="\t", index_col=0)
        mt = [g for g in pb.columns if g.startswith("MT-")]
        with tempfile.TemporaryDirectory() as td:
            pb.drop(columns=mt).to_csv(f"{td}/pb_noMT.tsv.gz", sep="\t")
            run(f"{RSCRIPT} {SCR}/01_limma_q1q2.R --pb {td}/pb_noMT.tsv.gz --det {DET} "
                f"--meta {DATA}/meta_T2000.tsv --thresh {THRESH} --out {td} "
                f"--tag _T2000_d25_noMT")
            for suffix in ("noMT", "noMTrenorm"):
                pd.read_csv(f"{td}/q1_raw_T2000_d25_noMT.tsv", sep="\t").to_csv(
                    f"{RES}/q1_raw_T2000_d25_{suffix}.tsv", sep="\t", index=False)
        print(f"[noMT] removed {len(mt)} mitochondrially encoded genes: {', '.join(mt)}")

    # --- summaries -----------------------------------------------------------------
    q1 = pd.read_csv(f"{RES}/q1_raw_T2000_d25.tsv", sep="\t").set_index("gene")
    nm = pd.read_csv(f"{RES}/q1_raw_T2000_d25_noMTrenorm.tsv", sep="\t").set_index("gene")

    def row(label, t):
        s = t[t["adj.P.Val"] < 0.05]
        return dict(run=label, n_tested=len(t), n_sig=len(s),
                    n_up_CS=int((s.logFC > 0).sum()),
                    median_abs_lfc=round(s.logFC.abs().median(), 3))

    pd.DataFrame([row("original", q1), row("MT-free, renormalized", nm)]).to_csv(
        f"{RES}/q1_noMT_summary.tsv", sep="\t", index=False)

    J = q1[["logFC", "adj.P.Val"]].join(
        nm[["logFC", "adj.P.Val"]], rsuffix="_noMT", how="inner").dropna()
    J.to_csv(f"{RES}/q1_noMT_comparison.tsv", sep="\t")

    # Restricted to main hits that are also testable in the comparator-only panel
    # contrast, so every row can be read against the directly measured panel effect.
    panel_testable = set(pd.read_csv(f"{RES}/ctrl_panelB_vs_panelA.tsv", sep="\t").gene)
    sig = [g for g in q1.index[q1["adj.P.Val"] < 0.05] if g in panel_testable]
    C = pd.DataFrame({"logFC": q1.loc[sig, "logFC"]})
    for name, tag in [("ctrl_panelB_vs_panelA", "pnl"), ("q1_panelmatched", "pm"),
                      ("q1_purityadj", "pu")]:
        t = pd.read_csv(f"{RES}/{name}.tsv", sep="\t").set_index("gene")
        C[f"logFC_{tag}"] = t.logFC.reindex(C.index)
        C[f"fdr_{tag}"] = t["adj.P.Val"].reindex(C.index)
    C = C.dropna()
    C.to_csv(f"{RES}/q1_control_comparison.tsv", sep="\t")

    print(pd.read_csv(f"{RES}/q1_noMT_summary.tsv", sep="\t").to_string(index=False))
    for tag, label in [("pm", "panel-matched"), ("pu", "purity-adjusted")]:
        same = int((np.sign(C.logFC) == np.sign(C[f"logFC_{tag}"])).sum())
        print(f"[control {label:<15}] {len(C)} main hits testable | "
              f"same direction {same} | still FDR<0.05 {int((C[f'fdr_{tag}'] < 0.05).sum())}")
    print(f"[panel effect] {int((C.fdr_pnl < 0.05).sum())} of {len(C)} main hits are "
          f"themselves differential between probe panels among comparators")


if __name__ == "__main__":
    main()
