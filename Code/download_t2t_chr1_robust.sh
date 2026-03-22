#!/usr/bin/env bash
# A more robust, single-file download script for the T2T chr1 VCF.
# It will try to use aria2c if available, otherwise fall back to a resilient curl.

set -euo pipefail

# --- Configuration ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_DIR="${SCRIPT_DIR}/../Data/T2T_CHM13"
BASE_URL="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/all_samples_3202"
VCF_FILE="1KGP.CHM13v2.0.chr1.recalibrated.snp_indel.pass.vcf.gz"
TBI_FILE="${VCF_FILE}.tbi"
LOG_FILE="${SCRIPT_DIR}/../Results/logs/robust_download.log"

mkdir -p "${OUTPUT_DIR}"
touch "${LOG_FILE}"

echo "--- Robust Download Script Started: $(date) ---" | tee -a "${LOG_FILE}"

download() {
    local url=$1
    local dest=$2
    local filename=$(basename "$dest")

    echo "Attempting to download ${filename}..." | tee -a "${LOG_FILE}"

    if [[ -f "$dest" ]]; then
        echo "${filename} already exists. Skipping." | tee -a "${LOG_FILE}"
        return 0
    fi

    if command -v aria2c &> /dev/null; then
        echo "Using aria2c for robust download." | tee -a "${LOG_FILE}"
        aria2c \
            --continue=true \
            --max-connection-per-server=8 \
            --split=8 \
            --min-split-size=1M \
            --retry-wait=5 \
            --max-tries=0 \
            -d "$(dirname "$dest")" \
            -o "$filename" \
            "$url" | tee -a "${LOG_FILE}"
    else
        echo "aria2c not found. Falling back to resilient curl." | tee -a "${LOG_FILE}"
        curl --fail --location --retry 20 --retry-delay 10 \
             -C - --connect-timeout 60 \
             -# -o "$dest" "$url" | tee -a "${LOG_FILE}"
    fi

    if [[ -f "$dest" ]]; then
        echo "SUCCESS: ${filename} downloaded." | tee -a "${LOG_FILE}"
        return 0
    else
        echo "ERROR: Failed to download ${filename}." | tee -a "${LOG_FILE}"
        return 1
    fi
}

# --- Main Execution ---
download "${BASE_URL}/${VCF_FILE}" "${OUTPUT_DIR}/${VCF_FILE}"
download "${BASE_URL}/${TBI_FILE}" "${OUTPUT_DIR}/${TBI_FILE}"

echo "--- Robust Download Script Finished: $(date) ---" | tee -a "${LOG_FILE}"
