#!/bin/bash
# Download the GRCh38 chrX VCF from 1000 Genomes EBI FTP
# Note: chrX uses eagle2-phased.v2 (not shapeit2-duohmm) — this is expected
# for sex chromosomes where duo-HMM pedigree phasing is not applicable.

set -euo pipefail

BASE_URL="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased"
OUTPUT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/Data/GRCh38"
VCF="CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.vcf.gz"
TBI="${VCF}.tbi"

mkdir -p "${OUTPUT_DIR}"

download_file() {
    local url="$1"
    local dest="$2"
    local tmp="${dest}.part"

    if [[ -s "${dest}" ]]; then
        echo "Already exists: $(basename ${dest}) ($(du -sh ${dest} | cut -f1))"
        return 0
    fi

    echo "Downloading $(basename ${dest})..."
    if curl -L -C - --progress-bar -o "${tmp}" "${url}"; then
        mv "${tmp}" "${dest}"
        echo "Done: $(basename ${dest}) ($(du -sh ${dest} | cut -f1))"
    else
        rm -f "${tmp}"
        echo "ERROR: failed to download $(basename ${dest})"
        exit 1
    fi
}

download_file "${BASE_URL}/${VCF}" "${OUTPUT_DIR}/${VCF}"
download_file "${BASE_URL}/${TBI}" "${OUTPUT_DIR}/${TBI}"

echo ""
echo "=== GRCh38 chrX download complete ==="
ls -lh "${OUTPUT_DIR}/${VCF}" "${OUTPUT_DIR}/${TBI}"
