#!/bin/bash
# Download 1000 Genomes Project VCFs aligned to GRCh38 (NYGC 30x high-coverage)
# Source: 1000 Genomes FTP / NYGC
# 3,202 samples, phased VCFs for all chromosomes (1-22, X)
# This script downloads files sequentially to avoid overwhelming the server.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_DIR="${SCRIPT_DIR}/../Data/GRCh38"
BASE_URL="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased"

mkdir -p "${OUTPUT_DIR}"

download_file() {
    local url="$1"
    local dest="$2"
    local tmp="${dest}.part"

    if [[ -s "${dest}" ]]; then
        echo "Skipping $(basename "${dest}") - already exists"
        return 0
    fi

    echo "Downloading $(basename "${dest}")..."
    if curl -L -C - -o "${tmp}" "${url}" && mv "${tmp}" "${dest}"; then
        echo "Successfully downloaded $(basename "${dest}")"
    else
        echo "Error downloading $(basename "${dest}")"
        rm -f "${tmp}"
        return 1
    fi
}

# Chromosomes 1-22 and X
CHROMS=({1..22} X)

for CHR in "${CHROMS[@]}"; do
    VCF_FILENAME="CCDG_14151_B01_GRM_WGS_2020-08-05_chr${CHR}.filtered.shapeit2-duohmm-phased.vcf.gz"
    TBI_FILENAME="${VCF_FILENAME}.tbi"
    
    VCF_URL="${BASE_URL}/${VCF_FILENAME}"
    TBI_URL="${BASE_URL}/${TBI_FILENAME}"
    
    VCF_DEST="${OUTPUT_DIR}/${VCF_FILENAME}"
    TBI_DEST="${OUTPUT_DIR}/${TBI_FILENAME}"

    download_file "${VCF_URL}" "${VCF_DEST}"
    download_file "${TBI_URL}" "${TBI_DEST}"
done

echo "All GRCh38 files downloaded successfully."
