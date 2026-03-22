#!/bin/bash
# Download 1000 Genomes Project VCFs aligned to GRCh38 (NYGC 30x high-coverage) from AWS S3
# Source: 1000 Genomes FTP / NYGC - Mirrored on AWS S3
# 3,202 samples, phased VCFs for all chromosomes (1-22, X)
# This script downloads files sequentially to avoid overwhelming the server.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_DIR="${SCRIPT_DIR}/../Data/GRCh38"
# Using the AWS S3 mirror for better reliability
BASE_URL="http://s3.amazonaws.com/1000genomes/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased"

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
    # Using aria2c for more robust downloading
    if aria2c -x 4 -s 4 --auto-file-renaming=false -c -d "$(dirname "${dest}")" -o "$(basename "${dest}")" "${url}"; then
        echo "Successfully downloaded $(basename "${dest}")"
    else
        echo "Error downloading $(basename "${dest}")"
        # aria2c creates a .aria2 file, remove it on failure
        rm -f "${dest}.aria2"
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
