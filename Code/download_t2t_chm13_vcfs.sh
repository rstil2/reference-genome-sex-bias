
#!/bin/bash
# Robust download of 1000 Genomes Project VCFs aligned to T2T-CHM13v2.0
# Source: T2T Consortium AWS bucket
# 3,202 samples, all chromosomes (1-22, X)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_DIR="${SCRIPT_DIR}/../Data/T2T_CHM13"
BASE_URL="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/all_samples_3202"
LOG_FILE="${OUTPUT_DIR}/download.log"
FAILED_FILE="${OUTPUT_DIR}/failed_downloads.txt"

mkdir -p "${OUTPUT_DIR}"
rm -f "$LOG_FILE" "$FAILED_FILE"

# Download a file robustly, with retries and logging
download_file() {
    local url="$1"
    local dest="$2"
    local tmp="${dest}.part"
    local max_retries=20
    local retry_count=0

    if [[ -s "${dest}" ]]; then
        echo "Skipping $(basename "${dest}") - already exists" | tee -a "$LOG_FILE"
        return 0
    fi

    while (( retry_count < max_retries )); do
        echo "Downloading $(basename "${dest}") (attempt $((retry_count+1))/$max_retries)..." | tee -a "$LOG_FILE"
        curl --fail --location --retry 10 --retry-delay 5 \
             -C - --connect-timeout 60 \
             -# -o "${tmp}" "${url}"

        if [[ -s "${tmp}" ]]; then
            mv "${tmp}" "${dest}"
            echo "SUCCESS: $(basename "${dest}")" | tee -a "$LOG_FILE"
            return 0
        else
            echo "WARNING: Downloaded file is empty: ${tmp}" | tee -a "$LOG_FILE"
            rm -f "${tmp}"
        fi
        ((retry_count++))
        sleep 2
    done
    echo "FAILED: $(basename "${dest}") after $max_retries attempts" | tee -a "$LOG_FILE"
    echo "${url}" >> "$FAILED_FILE"
    return 1
}


echo "Downloading T2T-CHM13 VCFs to ${OUTPUT_DIR}"
echo "Source: ${BASE_URL}"
echo "Log: ${LOG_FILE}"
echo "Failed downloads: ${FAILED_FILE}"
echo ""


# Download autosomes and X chromosome, skipping chr7 on first pass
for chr in {1..22} X; do
    if [[ "$chr" == "7" ]]; then
        continue
    fi
    VCF_FILE="1KGP.CHM13v2.0.chr${chr}.recalibrated.snp_indel.pass.vcf.gz"
    TBI_FILE="${VCF_FILE}.tbi"

    download_file "${BASE_URL}/${VCF_FILE}" "${OUTPUT_DIR}/${VCF_FILE}"
    download_file "${BASE_URL}/${TBI_FILE}" "${OUTPUT_DIR}/${TBI_FILE}"
done

# Now attempt chr7 at the end
chr="7"
VCF_FILE="1KGP.CHM13v2.0.chr${chr}.recalibrated.snp_indel.pass.vcf.gz"
TBI_FILE="${VCF_FILE}.tbi"
download_file "${BASE_URL}/${VCF_FILE}" "${OUTPUT_DIR}/${VCF_FILE}"
download_file "${BASE_URL}/${TBI_FILE}" "${OUTPUT_DIR}/${TBI_FILE}"

echo ""
echo "Download complete. Files saved to ${OUTPUT_DIR}"
echo "Total VCF files: $(ls -lh "${OUTPUT_DIR}"/*.vcf.gz 2>/dev/null | wc -l)"
echo "Total TBI files: $(ls -lh "${OUTPUT_DIR}"/*.tbi 2>/dev/null | wc -l)"
if [[ -s "$FAILED_FILE" ]]; then
    echo "Some files failed to download. See $FAILED_FILE for retry links."
else
    echo "All files downloaded successfully."
fi
