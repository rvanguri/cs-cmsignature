#!/usr/bin/env python3
"""Download open GEO datasets (Neyazi, Larson, Chin2022, Chin2021) to $RAW.

Reads an accessions TSV (cols: dataset<TAB>gse<TAB>supp_url_or_blank) and pulls the supplementary
matrices for each GSE. If supp_url is blank it queries GEO's FTP supplementary path by convention.

This is a thin, robust wrapper around `wget`; it does NOT hard-code per-sample URLs so it works as the
accession list evolves. Fill scripts/../raw/accessions.tsv before running 01_download.sbatch.

TODO(you): confirm each GSE's exact supplementary filename(s); GEO layouts vary by submitter.
"""
from __future__ import annotations
import argparse
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import log, ensure_dir, require  # noqa: E402

LOG = log("download_geo")


def geo_supp_ftp(gse: str) -> str:
    """Conventional GEO FTP supplementary directory for a GSE accession."""
    stub = gse[:-3] + "nnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{gse}/suppl/"


def fetch(url: str, dest_dir: str) -> None:
    ensure_dir(dest_dir)
    LOG.info("wget %s -> %s", url, dest_dir)
    # -r -np -nd: recurse one level into the suppl dir, no parent, flatten; -A restrict to data files
    cmd = ["wget", "-r", "-np", "-nd", "-e", "robots=off",
           "-A", "*.h5,*.tar,*.tar.gz,*.mtx.gz,*.tsv.gz,*.csv.gz,*.txt.gz",
           "-P", dest_dir, url]
    subprocess.run(cmd, check=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True,
                    help="accessions TSV: dataset<TAB>gse<TAB>supp_url(optional)")
    ap.add_argument("--out", required=True, help="raw/ root")
    args = ap.parse_args()
    require(args.manifest, "accessions TSV")

    with open(args.manifest) as fh:
        for ln in fh:
            ln = ln.strip()
            if not ln or ln.startswith("#"):
                continue
            parts = ln.split("\t")
            dataset, gse = parts[0], parts[1]
            url = parts[2] if len(parts) > 2 and parts[2] else geo_supp_ftp(gse)
            fetch(url, os.path.join(args.out, dataset))
    LOG.info("GEO downloads complete.")


if __name__ == "__main__":
    main()
