#!/usr/bin/env python3
"""Per-locus quantification QC for the genes tested in the spatial-first analysis.

Motivation. The comparator sections and the CS sections were generated on different
Visium chemistries in different labs. A gene can therefore differ between cohorts for
reasons that have nothing to do with biology: the probe set covers it differently, the
symbol is a readthrough locus whose reads are split between two annotations, or the
locus sits in a segmental duplication. This step flags each tested gene on five
independent criteria and marks a gene `qc_pass` only if it trips none of them.

Flags
  flag_source_specific_absence  detected (mean detection fraction >= 0.25) in >= 1 data
                                source and effectively absent (max < 0.02) in >= 1 other
  flag_length_deviation         cross-source expression deviation > 2 log2 units after
                                regressing out gene span (probe-panel length bias)
  flag_segdup                   gene span overlaps an hg38 segmental duplication
  flag_missing_annotation       symbol not resolvable to an Ensembl gene
  flag_readthrough              symbol is a hyphenated readthrough whose two components
                                are themselves genes in the reference symbol space

Inputs (all committed under data/spatial_first/)
  pseudobulk_T2000.tsv.gz            sections x genes raw pseudobulk counts
  meta_T2000.tsv                     section metadata (cohort, source, evaluable)
  detection_by_section_T2000.tsv.gz  per-section detection fraction
  ensembl_annotation_tested.tsv      cached Ensembl lookup for the tested genes
  segdups_hg38.tsv.gz                cached UCSC genomicSuperDups track
  reference_gene_symbols.txt         10x reference symbol space (readthrough detection)
Input from the DE step
  results/spatial_first/q2_raw_T2000_d25.tsv   defines the tested-gene universe

Pass --refresh to re-fetch the Ensembl and UCSC tables from the network instead of
using the committed caches; the caches exist so this step runs offline.

Output: results/spatial_first/q2_locus_qc.tsv
"""
import argparse
import json
import time

import numpy as np
import pandas as pd

DATA = "data/spatial_first"
RES = "results/spatial_first"
FLAGS = ["flag_source_specific_absence", "flag_length_deviation", "flag_segdup",
         "flag_missing_annotation", "flag_readthrough"]
SOURCES = ["foong", "new_visium", "kuppe"]


def fetch_ensembl(genes):
    import requests
    url = "https://rest.ensembl.org/lookup/symbol/homo_sapiens"
    info = {}
    for i in range(0, len(genes), 250):
        r = requests.post(url, headers={"Content-Type": "application/json",
                                        "Accept": "application/json"},
                          data=json.dumps({"symbols": genes[i:i + 250]}), timeout=180)
        r.raise_for_status()
        info.update(r.json())
        time.sleep(0.3)
    ann = pd.DataFrame([{"gene": g, "ens_id": v.get("id"), "chrom": v.get("seq_region_name"),
                         "start": v.get("start"), "end": v.get("end"), "strand": v.get("strand"),
                         "biotype": v.get("biotype"), "description": v.get("description") or ""}
                        for g, v in info.items() if isinstance(v, dict)])
    ann["span_kb"] = (ann.end - ann.start) / 1000
    return ann


def fetch_segdups():
    import requests
    rows = []
    for c in [str(i) for i in range(1, 23)] + ["X", "Y"]:
        r = requests.get("https://api.genome.ucsc.edu/getData/track",
                         params={"genome": "hg38", "track": "genomicSuperDups",
                                 "chrom": "chr" + c}, timeout=300)
        r.raise_for_status()
        rows += [{"chrom": c, "start": x["chromStart"], "end": x["chromEnd"],
                  "frac": x.get("fracMatch")} for x in r.json().get("genomicSuperDups", [])]
    return pd.DataFrame(rows)


def readthrough_map(refsym):
    """Symbols whose hyphenated readthrough partner is also in the reference space.

    MT- and HLA- prefixes are excluded: those hyphens are naming conventions, not
    readthrough annotations.
    """
    upper = set(pd.Index(refsym).str.upper())
    hyphenated = [s for s in upper
                  if "-" in s and not s.startswith("MT-") and not s.startswith("HLA-")]
    out = {}
    for s in hyphenated:
        parts = s.split("-")
        for k in range(1, len(parts)):
            a, b = "-".join(parts[:k]), "-".join(parts[k:])
            if a in upper and b in upper:
                out.setdefault(a, []).append(s)
                out.setdefault(b, []).append(s)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--refresh", action="store_true",
                    help="re-fetch Ensembl/UCSC instead of using the committed caches")
    ap.add_argument("--out", default=f"{RES}/q2_locus_qc.tsv")
    args = ap.parse_args()

    pb = pd.read_csv(f"{DATA}/pseudobulk_T2000.tsv.gz", sep="\t", index_col=0)
    meta = pd.read_csv(f"{DATA}/meta_T2000.tsv", sep="\t", index_col=0)
    dbs = pd.read_csv(f"{DATA}/detection_by_section_T2000.tsv.gz", sep="\t", index_col=0)
    refsym = [x.strip() for x in open(f"{DATA}/reference_gene_symbols.txt") if x.strip()]
    genes = pd.read_csv(f"{RES}/q2_raw_T2000_d25.tsv", sep="\t").gene.tolist()

    ev = meta[meta.evaluable.astype(str).isin(["True", "TRUE"])]

    if args.refresh:
        ann = fetch_ensembl(genes)
        segdup = fetch_segdups()
    else:
        ann = pd.read_csv(f"{DATA}/ensembl_annotation_tested.tsv", sep="\t")
        segdup = pd.read_csv(f"{DATA}/segdups_hg38.tsv.gz", sep="\t")

    # ---- source-specific absence -------------------------------------------------
    dbs = dbs.loc[[s for s in ev.index if s in dbs.index]]
    det = pd.DataFrame(index=pb.columns)
    for s in SOURCES:
        sub = dbs.loc[ev.source[ev.source == s].index]
        det[f"det_max_{s}"] = sub.max().reindex(det.index)
        det[f"det_mean_{s}"] = sub.mean().reindex(det.index)
    absent = det[[f"det_max_{s}" for s in SOURCES]] < 0.02
    present = det[[f"det_mean_{s}" for s in SOURCES]] >= 0.25
    det["n_sources_absent"] = absent.sum(1)
    det["n_sources_present"] = present.sum(1)
    det["flag_source_specific_absence"] = (det.n_sources_absent >= 1) & (det.n_sources_present >= 1)

    # ---- length-matched cross-source deviation -----------------------------------
    cpm = np.log2(pb.loc[ev.index].div(pb.loc[ev.index].sum(1), axis=0) * 1e6 + 1)
    srcmean = cpm.groupby(ev.source).mean().T
    A = ann.set_index("gene").reindex(srcmean.index)
    ok = A.span_kb.notna() & (srcmean.max(1) > 0)
    x = np.log10(A.span_kb[ok].values)
    cross = srcmean[ok].mean(1)
    dev = pd.DataFrame(index=srcmean.index[ok])
    for s in srcmean.columns:
        y = (srcmean.loc[ok, s] - cross).values
        dev[s] = y - np.polyval(np.polyfit(x, y, 1), x)
    dev["max_abs_dev"] = dev.abs().max(1)

    # ---- segmental-duplication overlap -------------------------------------------
    by_chrom = {str(c): g[["start", "end"]].values for c, g in segdup.groupby("chrom")}

    def overlaps(row):
        iv = by_chrom.get(str(row.chrom))
        if iv is None or pd.isna(row.start):
            return False
        return bool(((iv[:, 0] < row.end) & (iv[:, 1] > row.start)).any())

    A["flag_segdup"] = A.apply(overlaps, axis=1)

    # ---- assemble -----------------------------------------------------------------
    rt = readthrough_map(refsym)
    QC = pd.DataFrame(index=pd.Index(genes, name="gene")).join(
        A[["ens_id", "chrom", "start", "end", "biotype", "span_kb", "description", "flag_segdup"]])
    QC["flag_readthrough_partner"] = [";".join(rt.get(g, [])) for g in QC.index]
    QC["flag_readthrough"] = QC.flag_readthrough_partner != ""
    QC["max_abs_length_dev_log2"] = dev.max_abs_dev.reindex(QC.index)
    QC["flag_length_deviation"] = QC.max_abs_length_dev_log2 > 2
    QC = QC.join(det[["flag_source_specific_absence", "n_sources_absent"]
                     + [f"det_mean_{s}" for s in SOURCES]])
    QC["flag_missing_annotation"] = QC.ens_id.isna()
    QC[FLAGS] = QC[FLAGS].fillna(False).astype(bool)
    QC["n_flags"] = QC[FLAGS].sum(1)
    QC["qc_pass"] = QC.n_flags == 0
    QC.to_csv(args.out, sep="\t")
    print(f"[locus qc] {len(QC)} tested genes | qc_pass {int(QC.qc_pass.sum())}")
    print(QC[FLAGS].sum().to_string())


if __name__ == "__main__":
    main()
