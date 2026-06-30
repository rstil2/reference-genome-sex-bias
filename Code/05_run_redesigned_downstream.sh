#!/usr/bin/env bash
# 05_run_redesigned_downstream.sh
# Run the v2 downstream panel battery (reference-sensitive + negative controls).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
PROJECT_DIR="$(cd "${ROOT_DIR}/.." && pwd)"
VENV_PY="${PROJECT_DIR}/.venv/bin/python"
PY="${VENV_PY}"
export DYLD_LIBRARY_PATH="/opt/homebrew/lib:${DYLD_LIBRARY_PATH:-}"

TS="$(date -u +"%Y%m%dT%H%M%SZ")"
LOG="${ROOT_DIR}/Results/logs/redesigned_downstream_${TS}.log"
mkdir -p "${ROOT_DIR}/Results/logs"

log() { echo "[$(date -u +"%Y-%m-%dT%H:%M:%SZ")] $*" | tee -a "${LOG}"; }

log "=== Redesigned Downstream Panel Battery (v2) ==="
log "Study root: ${ROOT_DIR}"

"${PY}" "${SCRIPT_DIR}/10_downstream_panel_battery.py" \
    --study-root "${ROOT_DIR}" \
    --python "${PY}" \
    --seed 20260222 \
    --max-repeats 20 \
    2>&1 | tee -a "${LOG}"

log "=== Complete. Results in Results/metrics/downstream_v2/ ==="
log "Copy revision summaries to Revisions/Analysis/ as needed (Google Drive project root)."
