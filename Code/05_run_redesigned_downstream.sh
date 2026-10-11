#!/usr/bin/env bash
# 05_run_redesigned_downstream.sh
# Run the v2 downstream panel battery (reference-sensitive + negative controls).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
if [[ -x "${ROOT_DIR}/.venv/bin/python" ]]; then
  PY="${ROOT_DIR}/.venv/bin/python"
elif command -v python3 >/dev/null 2>&1; then
  PY="$(command -v python3)"
else
  echo "python3 not found" >&2
  exit 1
fi

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
