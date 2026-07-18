#!/usr/bin/env bash
# Download the Foong cardiac-sarcoidosis spatial (Visium) series from GEO into $FOONG_DIR.
# Accession is provided via $FOONG_GSE (e.g. GSE######). We do NOT hard-code it because the
# exact series could not be verified programmatically: set it from the paper's Data Availability.
#
# GEO FTP layout: https://ftp.ncbi.nlm.nih.gov/geo/series/<STUB>/<GSE>/suppl/
#   where <STUB> = GSE + (accession-number with the last 3 digits replaced by 'nnn')
# Visium series usually ship either a combined .h5ad, a *_RAW.tar of per-sample Space Ranger
# outputs (filtered_feature_bc_matrix.h5 + spatial/), or per-sample MTX triplets. We pull
# everything in suppl/ and untar any archives so the figures script can auto-detect the format.
set -euo pipefail
: "${FOONG_GSE:?set FOONG_GSE to the Foong Visium GEO accession, e.g. GSE211469}"
: "${FOONG_DIR:?set FOONG_DIR (destination dir)}"

acc="${FOONG_GSE}"
num="${acc#GSE}"
if (( ${#num} > 3 )); then stub="GSE${num:0:${#num}-3}nnn"; else stub="GSEnnn"; fi
url="https://ftp.ncbi.nlm.nih.gov/geo/series/${stub}/${acc}/suppl/"

mkdir -p "${FOONG_DIR}"
echo "[foong] downloading ${url} -> ${FOONG_DIR}"
# -r recursive, -np no parent, -nd flatten, -e robots=off, retry, continue partial
wget -r -np -nd -e robots=off --tries=3 --continue -P "${FOONG_DIR}" "${url}" \
  || { echo "[foong] wget failed: check FOONG_GSE=${acc} / network egress on this node"; exit 1; }

echo "[foong] extracting archives (RAW.tar, then any nested per-sample archives)"
shopt -s nullglob
# pass 1: top-level RAW tar(s)
for t in "${FOONG_DIR}"/*.tar "${FOONG_DIR}"/*.tar.gz "${FOONG_DIR}"/*.tgz; do
    echo "  untar $t"; tar -xf "$t" -C "${FOONG_DIR}" || echo "  (skip $t)"
done
# pass 2: any nested spatial tarballs unpacked by pass 1 (Space Ranger sometimes tars spatial/)
for t in "${FOONG_DIR}"/*_spatial.tar.gz "${FOONG_DIR}"/*spatial*.tar*; do
    echo "  untar nested $t"; tar -xf "$t" -C "${FOONG_DIR}" || echo "  (skip $t)"
done
echo "[foong] contents:"
ls -la "${FOONG_DIR}" | head -40 || true    # || true: head closes the pipe -> SIGPIPE would trip pipefail
echo "[foong] download complete"
