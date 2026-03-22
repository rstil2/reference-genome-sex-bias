#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
LOG_DIR="${ROOT_DIR}/Results/logs"
MANIFEST_DIR="${ROOT_DIR}/Results/manifests"

mkdir -p "${LOG_DIR}" "${MANIFEST_DIR}"

TS="$(date -u +"%Y%m%dT%H%M%SZ")"
OUT_TXT="${LOG_DIR}/status_snapshot_${TS}.txt"
OUT_JSON="${MANIFEST_DIR}/download_verification_${TS}.json"
LATEST_JSON="${MANIFEST_DIR}/download_verification_latest.json"

{
  echo "timestamp_utc=${TS}"
  echo "project_root=${ROOT_DIR}"
  echo "--- disk ---"
  df -h "${ROOT_DIR}" || true
  echo "--- data sizes ---"
  du -sh "${ROOT_DIR}/Data/T2T_CHM13" "${ROOT_DIR}/Data/GRCh38" "${ROOT_DIR}/Data/Metadata" || true
  echo "--- active partial files ---"
  ls -lh "${ROOT_DIR}/Data/T2T_CHM13"/*.part 2>/dev/null | head -n 3 || echo "no .part files in T2T_CHM13"
  ls -lh "${ROOT_DIR}/Data/GRCh38"/*.part 2>/dev/null | head -n 3 || echo "no .part files in GRCh38"
  echo "--- recent logs ---"
  ls -lt "${LOG_DIR}" | head -n 10 || true
} > "${OUT_TXT}"

/Users/stillwell/miniconda3/bin/conda run -p /Users/stillwell/miniconda3 --no-capture-output python \
  /Users/stillwell/.vscode/extensions/ms-python.python-2026.2.0-darwin-arm64/python_files/get_output_via_markers.py \
  "${ROOT_DIR}/Code/01_verify_downloads.py" \
  --t2t-dir "${ROOT_DIR}/Data/T2T_CHM13" \
  --grch38-dir "${ROOT_DIR}/Data/GRCh38" \
  --out-json "${OUT_JSON}" > /dev/null

cp "${OUT_JSON}" "${LATEST_JSON}"

echo "Status snapshot written: ${OUT_TXT}"
echo "Verification JSON written: ${OUT_JSON}"
