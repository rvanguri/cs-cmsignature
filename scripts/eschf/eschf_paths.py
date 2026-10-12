"""Input/output locations for the ESC HF cascade scripts.

ESCHF_PRIMARY=old (default): the original >=0.5 mm primary; reproduces the committed results/eschf tables.
ESCHF_PRIMARY=v2           : the all-cardiomyocyte-dominant-spot primary (handoff_figure_fixes s.0);
                             outputs go to results/eschf/primary_v2/ and fixed-count asserts are skipped.
"""
import os

MODE = os.environ.get("ESCHF_PRIMARY", "old")
assert MODE in ("old", "v2"), MODE
V2 = MODE == "v2"
OUT = "results/eschf/primary_v2" if V2 else "results/eschf"           # also holds q1_hits_annotated.tsv
RES = "results/eschf/primary_v2" if V2 else "results/spatial_first"   # q1 fit, LOO refits, snRNA per-gene
Q1 = f"{RES}/q1_raw_allcm.tsv" if V2 else f"{RES}/q1_raw_T2000_d25.tsv"
PROGRAM = f"{OUT}/geneset_primary_program_qcpass_down.txt" if V2 else "metadata/geneset_primary_program_qcpass_down.txt"
