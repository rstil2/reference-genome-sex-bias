#!/bin/bash
# Download T2T-CHM13 VCFs (chr1-22, X) using aria2c for fast, robust parallel downloads
# Requires: aria2c installed (brew install aria2 or apt-get install aria2)

OUTPUT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../Data/T2T_CHM13"
URL_LIST="${OUTPUT_DIR}/t2t_vcf_urls.txt"

mkdir -p "$OUTPUT_DIR"

BASE_URL="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/all_samples_3202"

# Generate URL list for chromosomes 1-22 and X
rm -f "$URL_LIST"
for chr in {1..22} X; do
    echo "$BASE_URL/1KGP.CHM13v2.0.chr${chr}.recalibrated.snp_indel.pass.vcf.gz" >> "$URL_LIST"
    echo "$BASE_URL/1KGP.CHM13v2.0.chr${chr}.recalibrated.snp_indel.pass.vcf.gz.tbi" >> "$URL_LIST"
done

# Download all files in parallel with aria2c
aria2c --input-file="$URL_LIST" \
       --dir="$OUTPUT_DIR" \
       --max-concurrent-downloads=8 \
       --continue=true \
       --retry-wait=5 \
       --max-tries=20 \
       --summary-interval=60 \
       --log="$OUTPUT_DIR/aria2_download.log" \
       --console-log-level=notice

echo "Download complete. Check $OUTPUT_DIR for files and aria2_download.log for details."
