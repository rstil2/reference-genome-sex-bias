#!/bin/bash
# Download 1000 Genomes Project sample metadata
# Contains: sample ID, sex, population, superpopulation, family info

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_DIR="${SCRIPT_DIR}/../Data/Metadata"
BASE_URL="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage"

mkdir -p "${OUTPUT_DIR}"

echo "Downloading sample metadata to ${OUTPUT_DIR}"
echo ""

# Main pedigree file with population info
PED_FILE="20130606_g1k_3202_samples_ped_population.txt"
echo "Downloading pedigree/population file..."
curl -# -o "${OUTPUT_DIR}/${PED_FILE}" "${BASE_URL}/${PED_FILE}"

# Also get the igsr population superpopulation mapping
SUPERPOP_URL="https://www.internationalgenome.org/api/beta/population/_search"
echo "Downloading population-superpopulation mapping..."
curl -s "https://www.internationalgenome.org/data-portal/api/population/_search?size=100" > "${OUTPUT_DIR}/population_superpopulation_mapping.json"

echo ""
echo "Download complete."
echo ""
echo "Files saved:"
ls -lh "${OUTPUT_DIR}/"

echo ""
echo "Sample metadata preview (first 5 lines):"
head -5 "${OUTPUT_DIR}/${PED_FILE}"
