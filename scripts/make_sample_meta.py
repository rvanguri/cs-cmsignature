#!/usr/bin/env python3
"""Generate a sample_meta.tsv skeleton from the staged raw/ folders.

Scans raw/<dataset>/<sample>/ and emits one row per sample, pre-filling the columns that are constant
per dataset (accession, chemistry, region, and disease/procurement where unambiguous) and INFERRING
disease/procurement from the sample folder name where the dataset is mixed (Chin2022 HCM/NF; Liu CS/ICM).
Anything it can't resolve is written as 'TODO' and listed in a REVIEW summary so you fill only those.

Reichart is intentionally skipped (its obs carries metadata; it isn't CellBendered per-sample).

Usage:
  python make_sample_meta.py --raw $RAW --out $RAW/sample_meta.tsv
Then open the file, fix every 'TODO', and sanity-check region/procurement for Neyazi (which varies).
"""
from __future__ import annotations
import argparse
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import log  # noqa: E402

LOG = log("make_sample_meta")

# constant per-dataset defaults; None => must be inferred from the sample name or left TODO
CFG = {
    "neyazi":   dict(accession="GSE319770", chemistry="3p", region="LV",  procurement="explant", disease="CS"),
    "larson":   dict(accession="GSE174691", chemistry="3p", region="IVS", procurement="myectomy", disease="HCM"),
    "chin2022": dict(accession="GSE181764", chemistry="3p", region="IVS", procurement=None,        disease=None),
    "chin2021": dict(accession="GSE161921", chemistry="3p", region="IVS", procurement="donor",    disease="NF"),
    "liu":      dict(accession="GSE205734", chemistry="5p", region="LV",  procurement="explant",  disease=None),
}
COLS = ["sample_id", "dataset", "accession", "chemistry", "procurement", "region", "disease", "sex", "individual"]


def infer(dataset: str, sample: str, cfg: dict) -> tuple[str, str, str]:
    """Return (disease, procurement, region), filling Nones from the sample-name tokens."""
    s = sample.lower()
    disease, proc, region = cfg["disease"], cfg["procurement"], cfg["region"]

    if dataset == "chin2022":               # mixed HCM (myectomy) / NF (donor)
        if any(t in s for t in ("hcm", "myectomy", "obstruct")):
            disease, proc = "HCM", "myectomy"
        elif any(t in s for t in ("nf", "donor", "control", "normal", "ctrl")):
            disease, proc = "NF", "donor"
        else:
            disease, proc = "TODO", "TODO"

    if dataset == "liu":                    # mixed CS / ICM
        if any(t in s for t in ("icm", "ischem", "ischaem")):
            disease = "ICM"
        elif any(t in s for t in ("cs", "sarc")):
            disease = "CS"
        else:
            disease = "TODO"

    if dataset == "neyazi":                 # CS, but procurement + region vary across samples
        if "vad" in s:
            proc = "VAD"
        elif any(t in s for t in ("autopsy", "postmortem", "pm")):
            proc = "autopsy"
        if any(t in s for t in ("_rv", "-rv", "rv_", "right")):
            region = "RV"
        elif any(t in s for t in ("_sp", "septum", "septal")):
            region = "SP"

    return disease or "TODO", proc or "TODO", region or "TODO"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    if os.path.exists(args.out):
        sys.exit(f"refusing to overwrite existing {args.out} (delete it first if you mean to regenerate)")

    rows, review = [], []
    for dataset, cfg in CFG.items():
        samples = sorted(os.path.basename(p.rstrip("/"))
                         for p in glob.glob(os.path.join(args.raw, dataset, "*/")))
        if not samples:
            LOG.warning("no sample folders under %s/%s: skipping", args.raw, dataset)
            continue
        for sample in samples:
            disease, proc, region = infer(dataset, sample, cfg)
            rows.append([sample, dataset, cfg["accession"], cfg["chemistry"],
                         proc, region, disease, "NA", sample])
            if "TODO" in (disease, proc):
                review.append(f"{dataset}/{sample}  disease={disease} procurement={proc}")
        LOG.info("%s: %d samples", dataset, len(samples))

    with open(args.out, "w") as f:
        f.write("# auto-generated skeleton: FIX every 'TODO' and verify Neyazi region/procurement, sex.\n")
        f.write("# Reichart intentionally absent (metadata in its own obs).\n")
        f.write("\t".join(COLS) + "\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")

    LOG.info("wrote %s (%d samples)", args.out, len(rows))
    if review:
        LOG.warning("---- %d rows need manual disease/procurement ----", len(review))
        for r in review:
            print("  REVIEW:", r)
    print("\nAlso confirm: Neyazi per-sample region (LV/RV/SP) + procurement (explant/VAD/autopsy), "
          "and 'sex' (M/F) wherever known.")


if __name__ == "__main__":
    main()
