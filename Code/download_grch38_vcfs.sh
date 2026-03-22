#!/bin/bash
# Download 1000 Genomes Project VCFs aligned to GRCh38 (NYGC 30x high-coverage)
# Source: 1000 Genomes FTP / NYGC
# 3,202 samples, phased VCFs for all chromosomes (1-22, X)

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
    # Use a robust curl command with resume and retry
    curl --fail --location --retry 10 --retry-delay 5 \
         -C - --connect-timeout 60 \
         -# -o "${tmp}" "${url}"

    if [[ ! -s "${tmp}" ]]; then
        echo "ERROR: Downloaded file is empty: ${tmp}"
        return 1
    fi

    mv "${tmp}" "${dest}"
}

echo "Downloading GRCh38 VCFs to ${OUTPUT_DIR}"
echo "Source: ${BASE_URL}"
echo ""

echo "Starting parallel download of all GRCh38 VCFs..."
echo "Log files for each download will be in ${OUTPUT_DIR}/logs"
mkdir -p "${OUTPUT_DIR}/logs"

# Download autosomes and X chromosome in parallel
for chr in {1..22} X; do
    VCF_FILE="CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.filtered.shapeit2-duohmm-phased.vcf.gz"
    TBI_FILE="${VCF_FILE}.tbi"
    
    (
        download_file "${BASE_URL}/${VCF_FILE}" "${OUTPUT_DIR}/${VCF_FILE}" &> "${OUTPUT_DIR}/logs/chr${chr}_vcf.log"
        download_file "${BASE_URL}/${TBI_FILE}" "${OUTPUT_DIR}/${TBI_FILE}" &> "${OUTPUT_DIR}/logs/chr${chr}_tbi.log"
    ) &
done

echo "Waiting for all parallel downloads to complete..."
wait


echo ""
echo "Download complete. Files saved to ${OUTPUT_DIR}"
echo "Total files:"
ls -lh "${OUTPUT_DIR}"/*.vcf.gz 2>/dev/null | wc -l
