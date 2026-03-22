#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

META_CSV="${ROOT_DIR}/Results/manifests/clean_metadata_v1.csv"
SPLIT_CSV="${ROOT_DIR}/Results/splits/split_manifest_v1.csv"

G_VCF="${ROOT_DIR}/Data/GRCh38/CCDG_14151_B01_GRM_WGS_2020-08-05_chr1.filtered.shapeit2-duohmm-phased.vcf.gz"
T_VCF="${ROOT_DIR}/Data/T2T_CHM13/1KGP.CHM13v2.0.chr1.recalibrated.snp_indel.pass.vcf.gz"

OUT_DIR="${ROOT_DIR}/Results/metrics/chr1_pilot"
PRED_DIR="${ROOT_DIR}/Results/predictions/chr1_pilot"
MANIFEST_DIR="${ROOT_DIR}/Results/manifests"
mkdir -p "${OUT_DIR}" "${PRED_DIR}" "${MANIFEST_DIR}"

if [[ ! -f "${G_VCF}" || ! -f "${T_VCF}" ]]; then
  echo "ERROR: chr1 VCF files not ready yet."
  echo "Missing GRCh38? $([[ -f "${G_VCF}" ]] && echo no || echo yes)"
  echo "Missing T2T? $([[ -f "${T_VCF}" ]] && echo no || echo yes)"
  exit 1
fi

PY="/Users/stillwell/projects/Project 33 - Bias in Reference Genomes/.venv/bin/python"
# Ensure the homebrew libomp (with required symbols) is found before the older /usr/local one
export DYLD_LIBRARY_PATH="/opt/homebrew/lib:${DYLD_LIBRARY_PATH:-}"

"${PY}" "${ROOT_DIR}/Code/02b_extract_variant_features.py" \
  --vcf "${G_VCF}" \
  --metadata "${META_CSV}" \
  --reference-label GRCh38 \
  --max-variants 2000 \
  --maf-threshold 0.01 \
  --call-rate-threshold 0.95 \
  --out-features-csv "${MANIFEST_DIR}/grch38_chr1_features_raw.csv" \
  --out-summary-json "${MANIFEST_DIR}/grch38_chr1_features_summary.json"

"${PY}" "${ROOT_DIR}/Code/02b_extract_variant_features.py" \
  --vcf "${T_VCF}" \
  --metadata "${META_CSV}" \
  --reference-label T2T_CHM13 \
  --max-variants 2000 \
  --maf-threshold 0.01 \
  --call-rate-threshold 0.95 \
  --out-features-csv "${MANIFEST_DIR}/t2t_chr1_features_raw.csv" \
  --out-summary-json "${MANIFEST_DIR}/t2t_chr1_features_summary.json"

"${PY}" "${ROOT_DIR}/Code/02c_align_reference_features.py" \
  --grch38-features "${MANIFEST_DIR}/grch38_chr1_features_raw.csv" \
  --t2t-features "${MANIFEST_DIR}/t2t_chr1_features_raw.csv" \
  --min-common-patients 100 \
  --min-features-per-ref 20 \
  --out-grch38-aligned "${MANIFEST_DIR}/grch38_chr1_features_aligned.csv" \
  --out-t2t-aligned "${MANIFEST_DIR}/t2t_chr1_features_aligned.csv" \
  --out-summary-json "${MANIFEST_DIR}/chr1_alignment_summary.json"

"${PY}" "${ROOT_DIR}/Code/03_train_models.py" \
  --features-csv "${MANIFEST_DIR}/grch38_chr1_features_aligned.csv" \
  --split-manifest "${SPLIT_CSV}" \
  --model logistic_regression \
  --max-repeats 10 \
  --metrics-out "${OUT_DIR}/grch38_lr_metrics.csv" \
  --predictions-out "${PRED_DIR}/grch38_lr_predictions.csv" \
  --summary-out "${OUT_DIR}/grch38_lr_summary.json"

"${PY}" "${ROOT_DIR}/Code/03_train_models.py" \
  --features-csv "${MANIFEST_DIR}/t2t_chr1_features_aligned.csv" \
  --split-manifest "${SPLIT_CSV}" \
  --model logistic_regression \
  --max-repeats 10 \
  --metrics-out "${OUT_DIR}/t2t_lr_metrics.csv" \
  --predictions-out "${PRED_DIR}/t2t_lr_predictions.csv" \
  --summary-out "${OUT_DIR}/t2t_lr_summary.json"

"${PY}" "${ROOT_DIR}/Code/04_analyze_results.py" \
  --grch38-metrics "${OUT_DIR}/grch38_lr_metrics.csv" \
  --t2t-metrics "${OUT_DIR}/t2t_lr_metrics.csv" \
  --bootstrap-iters 2000 \
  --permutation-iters 2000 \
  --out-csv "${OUT_DIR}/delta_delta_by_repeat.csv" \
  --out-json "${OUT_DIR}/delta_delta_summary.json"

echo "chr1 pilot complete. Outputs in: ${OUT_DIR} and ${PRED_DIR}"
