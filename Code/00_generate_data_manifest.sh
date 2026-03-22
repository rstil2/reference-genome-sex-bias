#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
DATA_DIR="${ROOT_DIR}/Data"
OUT_DIR="${ROOT_DIR}/Results/manifests"

mkdir -p "${OUT_DIR}"

TIMESTAMP="$(date -u +"%Y%m%dT%H%M%SZ")"
MANIFEST_TSV="${OUT_DIR}/data_manifest_${TIMESTAMP}.tsv"
MANIFEST_LATEST="${OUT_DIR}/data_manifest_latest.tsv"

{
  printf "path\tsha256\tsize_bytes\tmodified_utc\n"

  find "${DATA_DIR}" -type f \( \
    -name "*.vcf" -o -name "*.vcf.gz" -o -name "*.tbi" -o -name "*.tsv" -o -name "*.txt" -o -name "*.csv" \
  \) | sort | while IFS= read -r file; do
    rel_path="${file#${ROOT_DIR}/}"
    hash="$(shasum -a 256 "${file}" | awk '{print $1}')"
    size="$(stat -f "%z" "${file}")"
    mtime="$(date -u -r "$(stat -f "%m" "${file}")" +"%Y-%m-%dT%H:%M:%SZ")"
    printf "%s\t%s\t%s\t%s\n" "${rel_path}" "${hash}" "${size}" "${mtime}"
  done
} > "${MANIFEST_TSV}"

cp "${MANIFEST_TSV}" "${MANIFEST_LATEST}"

echo "Manifest written: ${MANIFEST_TSV}"
echo "Latest symlink copy: ${MANIFEST_LATEST}"
