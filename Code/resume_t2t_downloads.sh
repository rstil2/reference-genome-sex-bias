#!/bin/bash
# Resume incomplete T2T-CHM13 VCF downloads using aria2c.
#
# Context: A previous aria2c session downloaded chr1-12 completely but was
# interrupted for chr13-22 and chrX. The .aria2 control files are still present
# alongside the partial .vcf.gz files, so aria2c will resume from where it left off.
#
# Files skipped automatically (already complete, no .aria2 control file):
#   chr1-12 .vcf.gz and all .tbi files
#
# Usage: bash resume_t2t_downloads.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_DIR="${SCRIPT_DIR}/../Data/T2T_CHM13"
URL_LIST="${OUTPUT_DIR}/t2t_vcf_urls.txt"
BASE_URL="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/all_samples_3202"

# Timestamped log file to avoid appending to the existing 47 GB aria2_download.log
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="${OUTPUT_DIR}/aria2_resume_${TIMESTAMP}.log"

mkdir -p "$OUTPUT_DIR"

# Regenerate the URL list (chr1-22, X — both .vcf.gz and .tbi)
rm -f "$URL_LIST"
for chr in {1..22} X; do
    echo "$BASE_URL/1KGP.CHM13v2.0.chr${chr}.recalibrated.snp_indel.pass.vcf.gz" >> "$URL_LIST"
    echo "$BASE_URL/1KGP.CHM13v2.0.chr${chr}.recalibrated.snp_indel.pass.vcf.gz.tbi" >> "$URL_LIST"
done

echo "Resuming downloads to: $OUTPUT_DIR"
echo "Log: $LOG_FILE"
echo ""
echo "Files needing completion (have .aria2 control files):"
ls "$OUTPUT_DIR"/*.aria2 2>/dev/null | sed 's/.*\//  /' || echo "  (none — all downloads may already be complete)"
echo ""

# Resume all downloads.
# --continue=true  : resume partial downloads using existing .aria2 control files
# --auto-file-renaming=false : don't rename files if they already exist
# Files that are already complete (no .aria2 control file) are verified and skipped.
aria2c \
    --input-file="$URL_LIST" \
    --dir="$OUTPUT_DIR" \
    --max-concurrent-downloads=4 \
    --continue=true \
    --auto-file-renaming=false \
    --retry-wait=10 \
    --max-tries=30 \
    --summary-interval=120 \
    --log="$LOG_FILE" \
    --log-level=notice \
    --console-log-level=notice

echo ""
echo "Done. Checking final state..."
echo "Complete VCF files:  $(ls "$OUTPUT_DIR"/*.vcf.gz 2>/dev/null | grep -v '\.part$' | wc -l | tr -d ' ')"
echo "Complete TBI files:  $(ls "$OUTPUT_DIR"/*.tbi 2>/dev/null | wc -l | tr -d ' ')"
REMAINING=$(ls "$OUTPUT_DIR"/*.aria2 2>/dev/null | wc -l | tr -d ' ')
if [[ "$REMAINING" -gt 0 ]]; then
    echo "Still incomplete (.aria2 files remain): $REMAINING"
    ls "$OUTPUT_DIR"/*.aria2 | sed 's/.*\//  /'
else
    echo "All downloads complete — no .aria2 control files remaining."
fi
