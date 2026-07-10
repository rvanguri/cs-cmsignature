#!/usr/bin/env python3
"""Fetch the HeartMap harmonized atlas from Broad Single Cell Portal (SCP3689).

SCP downloads usually require an authenticated `curl`/bulk-download token (SCP "Bulk download" gives a
ready-made curl command with an auth_code). This script can't carry your credentials, so by default it
prints instructions and exits 0 (non-fatal) — the Plan's default arm mode is atlas-transfer, but the
pipeline still runs on open arms if the atlas isn't present yet.

If you have an SCP bulk-download curl/auth, pass it via --scp_curl_file (a file containing the command)
and this will execute it into --out.

TODO(you): obtain the SCP3689 bulk-download command (Broad SCP -> Download -> Bulk download).
"""
from __future__ import annotations
import argparse
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import log, ensure_dir  # noqa: E402

LOG = log("download_heartmap")
INSTRUCTIONS = """\
SCP3689 (HeartMap) not auto-downloaded. To get it:
  1) Visit https://singlecell.broadinstitute.org/single_cell/study/SCP3689
  2) Download -> Bulk download -> copy the generated `curl ... auth_code=...` command
  3) Save that command to a file and re-run:
       download_heartmap_scp3689.py --out <dir> --scp_curl_file <file>
  (Atlas-transfer is the Plan default; controlled raw via dbGaP phs001539 only if you flip an arm.)
"""


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--scp_curl_file", default="",
                    help="file containing the SCP bulk-download curl command (with auth_code)")
    args = ap.parse_args()
    ensure_dir(args.out)

    if not args.scp_curl_file:
        LOG.warning(INSTRUCTIONS)
        sys.exit(0)  # non-fatal: open arms still run

    with open(args.scp_curl_file) as fh:
        cmd = fh.read().strip()
    LOG.info("Running SCP bulk download into %s", args.out)
    subprocess.run(cmd, shell=True, check=True, cwd=args.out)
    LOG.info("HeartMap atlas download complete.")


if __name__ == "__main__":
    main()
