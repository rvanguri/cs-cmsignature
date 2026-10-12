#!/usr/bin/env python3
"""Regenerate derived input files that have no writer elsewhere in the repo.

Each rule was checked against the committed files: `check` reproduces every file exactly.

  data/spatial_first/meta_T2000.tsv          region_selection.tsv rows with T == 2000
  data/spatial_first/meta_panel.tsv          evaluable rows of meta_T2000 + probe panel;
                                             group = panel for comparators
  data/spatial_first/meta_q1_panelB.tsv      group = CS for CS, COMPB for panel-B comparators
  data/spatial_first/meta_q1_panelmatched.tsv group = CS for CS, COMPA for panel-A comparators
  metadata/tested_genes_985.txt              genes tested in q1_raw_T2000_d25.tsv (sensitivity arm)
  metadata/geneset_primary_program_qcpass_down.txt
                                             results/eschf/q1_hits_annotated.tsv rows with
                                             qc_pass and direction == "down"

Probe panel is a fixed per-section attribute (Foong CS: panelA; comparators as listed below);
it does not depend on the data.

Usage (repo root):
  python scripts/eschf/derived_inputs.py meta     # after the discovery export
  python scripts/eschf/derived_inputs.py lists    # after 00_inputs.py
  python scripts/eschf/derived_inputs.py check    # compare rules against files on disk
"""
import os
import sys

import pandas as pd

D = "data/spatial_first"
PANEL = {"foong": "panelA",
         "ARVC_1": "panelB", "ARVC_2": "panelB", "DC_LMNA": "panelB", "LMNA": "panelB",
         "Chagas": "panelA", "HCM": "panelA"}
META_COLS = ["section", "cohort", "disease", "source", "n_selected", "evaluable"]


def panel_of(r):
    return PANEL[r.section] if r.cohort == "COMPARATOR" else PANEL[r.source]


def build_meta():
    rs = pd.read_csv(f"{D}/region_selection.tsv", sep="\t")
    m = rs[rs["T"] == 2000][META_COLS].reset_index(drop=True)
    out = {"meta_T2000.tsv": m}
    e = m[m.evaluable.astype(str) == "True"].reset_index(drop=True).copy()
    e["panel"] = [panel_of(r) for r in e.itertuples()]
    comp = e.cohort == "COMPARATOR"
    p = e.copy(); p["group"] = p.panel.where(comp)
    b = e.copy(); b["group"] = None
    b.loc[b.cohort == "CS", "group"] = "CS"; b.loc[comp & (b.panel == "panelB"), "group"] = "COMPB"
    a = e.copy(); a["group"] = None
    a.loc[a.cohort == "CS", "group"] = "CS"; a.loc[comp & (a.panel == "panelA"), "group"] = "COMPA"
    out.update({"meta_panel.tsv": p, "meta_q1_panelB.tsv": b, "meta_q1_panelmatched.tsv": a})
    return out


def build_lists():
    q1 = pd.read_csv("results/spatial_first/q1_raw_T2000_d25.tsv", sep="\t")
    h = pd.read_csv("results/eschf/q1_hits_annotated.tsv", sep="\t")
    prog = h[(h.qc_pass == True) & (h.direction == "down")].gene.tolist()  # noqa: E712
    return {"metadata/tested_genes_985.txt": q1.gene.tolist(),
            "metadata/geneset_primary_program_qcpass_down.txt": prog}


def read_list(p):
    return [l.strip() for l in open(p) if l.strip() and not l.startswith("#")]


def main(mode):
    if mode == "check":
        for f, df in build_meta().items():
            old = pd.read_csv(f"{D}/{f}", sep="\t")
            new = df.copy(); new.to_csv("/tmp/_v3chk.tsv", sep="\t", index=False)
            new = pd.read_csv("/tmp/_v3chk.tsv", sep="\t")
            print(f"{f}: identical={old.equals(new)} rows={len(old)}")
        for f, g in build_lists().items():
            print(f"{f}: identical={read_list(f) == g} n={len(g)}")
        return
    if mode == "meta":
        for f, df in build_meta().items():
            df.to_csv(f"{D}/{f}", sep="\t", index=False)
            print(f"[meta] {f}: {len(df)} rows")
    elif mode == "lists":
        L = build_lists()
        for f, g in L.items():
            old_hdr = open(f).readline() if os.path.exists(f) else ""
            with open(f, "w") as fh:
                if f.endswith("qcpass_down.txt"):
                    fh.write("# QC-passing hits lower in CS (FDR<5%, q1_raw_T2000_d25 joined to "
                             f"q2_locus_qc); {len(g)} genes\n")
                elif old_hdr.startswith("#"):
                    fh.write(f"# genes tested in q1_raw_T2000_d25.tsv; {len(g)} genes\n")
                fh.write("\n".join(g) + "\n")
            print(f"[lists] {f}: {len(g)}")
    else:
        sys.exit(__doc__)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "")
