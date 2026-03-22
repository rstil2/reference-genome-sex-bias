#!/bin/bash
# Download 1000 Genomes Project VCFs aligned to T2T-CHM13v2.0 sequentially
# Source: T2T Consortium AWS bucket
# 3,202 samples, all chromosomes (1-22, X)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_DIR="${SCRIPT_DIR}/../Data/T2T_CHM13"
BASE_URL="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/all_samples_3202"

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
    if aria2c -x 4 -s 4 --auto-file-renaming=false -c -d "$(dirname "${dest}")" -o "$(basename "${dest}")" "${url}"; then
        echo "Successfully downloaded $(basename "${dest}")"
    else
        echo "Error downloading $(basename "${dest}")"
        rm -f "${dest}.aria2"
        return 1
    fi
}

echo "Downloading T2T-CHM13 VCFs to ${OUTPUT_DIR}"
echo "Source: ${BASE_URL}"
echo ""

# Download autosomes and X chromosome
for chr in {1..22} X; do
    VCF_FILE="1KGP.CHM13v2.0.chr${chr}.recalibrated.snp_indel.pass.vcf.gz"
    TBI_FILE="${VCF_FILE}.tbi"
    
    download_file "${BASE_URL}/${VCF_FILE}" "${OUTPUT_DIR}/${VCF_FILE}"
    download_file "${BASE_URL}/${TBI_FILE}" "${OUTPUT_DIR}/${TBI_FILE}"
done

echo "All T2T-CHM13 files downloaded successfully."
