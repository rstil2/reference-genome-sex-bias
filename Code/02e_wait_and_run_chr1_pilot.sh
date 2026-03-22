#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
LOG_DIR="${ROOT_DIR}/Results/logs"
mkdir -p "${LOG_DIR}"

G_VCF="${ROOT_DIR}/Data/GRCh38/CCDG_14151_B01_GRM_WGS_2020-08-05_chr1.filtered.shapeit2-duohmm-phased.vcf.gz"
T_VCF="${ROOT_DIR}/Data/T2T_CHM13/1KGP.CHM13v2.0.chr1.recalibrated.snp_indel.pass.vcf.gz"

TS="$(date -u +"%Y%m%dT%H%M%SZ")"
WATCH_LOG="${LOG_DIR}/chr1_watch_and_run_${TS}.log"

{
  echo "[$(date -u +"%Y-%m-%dT%H:%M:%SZ")] Watcher started"
  echo "Waiting for:"
  echo "  - ${G_VCF}"
  echo "  - ${T_VCF}"

  while true; do
    g_ready=0
    t_ready=0

    [[ -f "${G_VCF}" ]] && g_ready=1
    [[ -f "${T_VCF}" ]] && t_ready=1

    echo "[$(date -u +"%Y-%m-%dT%H:%M:%SZ")] status: GRCh38=${g_ready} T2T=${t_ready}"

    if [[ ${g_ready} -eq 1 && ${t_ready} -eq 1 ]]; then
      echo "[$(date -u +"%Y-%m-%dT%H:%M:%SZ")] Both chr1 VCFs ready. Running pilot pipeline."
      "${ROOT_DIR}/Code/02d_run_chr1_pilot_pipeline.sh"
      echo "[$(date -u +"%Y-%m-%dT%H:%M:%SZ")] Pilot pipeline completed."
      break
    fi

    sleep 300
  done
} >> "${WATCH_LOG}" 2>&1

echo "Watcher log: ${WATCH_LOG}"
