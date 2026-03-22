#!/bin/bash

# Set the base directory for data
DATA_DIR="../Data/GRCh38"
mkdir -p "$DATA_DIR"

# Base URL for the 1000 Genomes Project data on AWS S3
# This path is derived from the FTP path: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased_SNV_INDEL_SV/
BASE_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased_SNV_INDEL_SV"

echo "Downloading GRCh38 VCFs to $DATA_DIR"
echo "Source: $BASE_URL"

# Chromosomes to download
CHROMOSOMES=$(seq 1 22)
CHROMOSOMES+=" X"

# Loop through chromosomes and download each file
for CHR in $CHROMOSOMES; do
    # Construct the full URL for the file
    FILE_NAME="CCDG_14151_B01_GRM_WGS_2020-08-05_chr${CHR}.filtered.shapeit2-duohmm-phased.vcf.gz"
    REMOTE_URL="$BASE_URL/$FILE_NAME"
    LOCAL_FILE="$DATA_DIR/$FILE_NAME"
    
    echo "Downloading $FILE_NAME..."
    
    # Use aria2c for robust downloading
    # -c: continue getting a partially downloaded file
    # -x 4: use 4 connections per server
    # -s 4: split download to 4 files
    # --auto-file-renaming=false: do not rename file if it already exists
    aria2c -c -x 4 -s 4 --auto-file-renaming=false -o "$LOCAL_FILE" "$REMOTE_URL"
    
    if [ $? -ne 0 ]; then
        echo "Error downloading $FILE_NAME. Please check the URL and your connection."
        exit 1
    fi
done

echo "All GRCh38 VCF files downloaded successfully."
