#!/usr/bin/env bash
# 05_run_full_genome_pipeline.sh
# -----------------------------------------------------------------------
# Full Project 33 pipeline: generates splits, extracts genome-wide
# variant features from both references, aligns, trains all 5 model
# families, and computes the primary ΔΔ estimand with bootstrap CIs.
#
# Usage:
#   bash Study_v2_Real_Data/Code/05_run_full_genome_pipeline.sh [OPTIONS]
#
# Options (all optional, defaults shown):
#   --variants-per-chr N   Max variants extracted per chromosome (default: 500)
#   --maf FLOAT            MAF filter threshold (default: 0.01)
#   --call-rate FLOAT      Call-rate filter threshold (default: 0.95)
#   --repeats N            Number of split repeats (default: 20)
#   --seed N               Master random seed (default: 20260222)
#   --models "m1 m2 ..."   Space-separated model list (default: all 5)
#   --skip-phase1          Skip metadata/split generation (reuse existing)
#   --skip-phase2          Skip feature extraction (reuse existing)
#   --chr-list "c1 c2..."  Override chromosome list (default: 1-22 X)
# -----------------------------------------------------------------------

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
CODE_DIR="${SCRIPT_DIR}"
DATA_DIR="${ROOT_DIR}/Data"
RESULTS_DIR="${ROOT_DIR}/Results"
MANIFESTS_DIR="${RESULTS_DIR}/manifests"
SPLITS_DIR="${RESULTS_DIR}/splits"
METRICS_DIR="${RESULTS_DIR}/metrics"
PREDS_DIR="${RESULTS_DIR}/predictions"
CONFIGS_DIR="${RESULTS_DIR}/configs"
LOGS_DIR="${RESULTS_DIR}/logs"
FEATURES_DIR="${RESULTS_DIR}/features"

mkdir -p "${MANIFESTS_DIR}" "${SPLITS_DIR}" "${METRICS_DIR}" "${PREDS_DIR}" \
         "${CONFIGS_DIR}" "${LOGS_DIR}" "${FEATURES_DIR}"

# ---- Defaults ----
VARIANTS_PER_CHR=500
MAF=0.01
CALL_RATE=0.95
REPEATS=20
SEED=20260222
MODELS="logistic_regression random_forest xgboost svm mlp"
SKIP_PHASE1=0
SKIP_PHASE2=0
CHR_LIST="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X"

# ---- Argument parsing ----
while [[ $# -gt 0 ]]; do
    case "$1" in
        --variants-per-chr) VARIANTS_PER_CHR="$2"; shift 2 ;;
        --maf) MAF="$2"; shift 2 ;;
        --call-rate) CALL_RATE="$2"; shift 2 ;;
        --repeats) REPEATS="$2"; shift 2 ;;
        --seed) SEED="$2"; shift 2 ;;
        --models) MODELS="$2"; shift 2 ;;
        --skip-phase1) SKIP_PHASE1=1; shift ;;
        --skip-phase2) SKIP_PHASE2=1; shift ;;
        --chr-list) CHR_LIST="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

VENV_PY="/Users/stillwell/projects/Project 33 - Bias in Reference Genomes/.venv/bin/python"
PY="${VENV_PY}"
# Ensure the homebrew libomp (with required symbols) is found before the older /usr/local one
export DYLD_LIBRARY_PATH="/opt/homebrew/lib:${DYLD_LIBRARY_PATH:-}"
TS="$(date -u +"%Y%m%dT%H%M%SZ")"
RUN_LOG="${LOGS_DIR}/full_pipeline_${TS}.log"

log() { echo "[$(date -u +"%Y-%m-%dT%H:%M:%SZ")] $*" | tee -a "${RUN_LOG}"; }

log "=== Project 33 Full Genome Pipeline ==="
log "Root:              ${ROOT_DIR}"
log "Variants/chr:      ${VARIANTS_PER_CHR}"
log "MAF threshold:     ${MAF}"
log "Call-rate:         ${CALL_RATE}"
log "Repeats:           ${REPEATS}"
log "Seed:              ${SEED}"
log "Models:            ${MODELS}"
log "Chromosomes:       ${CHR_LIST}"

# Save config
cat > "${CONFIGS_DIR}/run_config_${TS}.json" <<EOF
{
  "generated_utc": "${TS}",
  "variants_per_chr": ${VARIANTS_PER_CHR},
  "maf_threshold": ${MAF},
  "call_rate_threshold": ${CALL_RATE},
  "repeats": ${REPEATS},
  "seed": ${SEED},
  "models": "$(echo ${MODELS} | tr ' ' ',')",
  "chromosomes": "$(echo ${CHR_LIST} | tr ' ' ',')"
}
EOF

# -----------------------------------------------------------------------
# PHASE 1: Metadata normalization and split manifest generation
# -----------------------------------------------------------------------
META_RAW="${DATA_DIR}/Metadata/20130606_g1k_3202_samples_ped_population.txt"
META_CLEAN="${MANIFESTS_DIR}/clean_metadata_v1.csv"
SPLIT_CSV="${SPLITS_DIR}/split_manifest_v1.csv"
PREPROCESS_SUMMARY="${MANIFESTS_DIR}/preprocess_summary_v1.json"

if [[ ${SKIP_PHASE1} -eq 0 ]]; then
    log "--- Phase 1: Generating clean metadata and split manifest ---"
    "${PY}" "${CODE_DIR}/02_preprocess.py" \
        --metadata "${META_RAW}" \
        --clean-metadata-out "${META_CLEAN}" \
        --split-manifest-out "${SPLIT_CSV}" \
        --summary-out "${PREPROCESS_SUMMARY}" \
        --train-frac 0.70 \
        --val-frac 0.15 \
        --test-frac 0.15 \
        --seed "${SEED}" \
        --repeats "${REPEATS}" \
        | tee -a "${RUN_LOG}"
    log "Phase 1 complete."
else
    log "--- Phase 1 skipped (--skip-phase1) ---"
    if [[ ! -f "${META_CLEAN}" || ! -f "${SPLIT_CSV}" ]]; then
        log "ERROR: --skip-phase1 set but required files are missing:"
        log "  ${META_CLEAN}"
        log "  ${SPLIT_CSV}"
        exit 1
    fi
fi

# -----------------------------------------------------------------------
# PHASE 2: Per-chromosome feature extraction + genome-wide merge
# -----------------------------------------------------------------------
GRCH38_MERGED="${FEATURES_DIR}/grch38_genome_features_raw.csv"
T2T_MERGED="${FEATURES_DIR}/t2t_genome_features_raw.csv"
GRCH38_ALIGNED="${FEATURES_DIR}/grch38_genome_features_aligned.csv"
T2T_ALIGNED="${FEATURES_DIR}/t2t_genome_features_aligned.csv"
ALIGN_SUMMARY="${MANIFESTS_DIR}/genome_alignment_summary.json"

if [[ ${SKIP_PHASE2} -eq 0 ]]; then
    log "--- Phase 2: Per-chromosome variant feature extraction ---"
    G_CHR_CSVS=()
    T_CHR_CSVS=()

    for CHR in ${CHR_LIST}; do
        # GRCh38 chrX uses eagle2-phased.v2 naming; all others use shapeit2-duohmm-phased
        if [[ "${CHR}" == "X" ]]; then
            G_VCF="${DATA_DIR}/GRCh38/CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.vcf.gz"
        else
            G_VCF="${DATA_DIR}/GRCh38/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${CHR}.filtered.shapeit2-duohmm-phased.vcf.gz"
        fi
        T_VCF="${DATA_DIR}/T2T_CHM13/1KGP.CHM13v2.0.chr${CHR}.recalibrated.snp_indel.pass.vcf.gz"

        G_OUT="${FEATURES_DIR}/grch38_chr${CHR}_features.csv"
        T_OUT="${FEATURES_DIR}/t2t_chr${CHR}_features.csv"
        G_SUM="${MANIFESTS_DIR}/grch38_chr${CHR}_features_summary.json"
        T_SUM="${MANIFESTS_DIR}/t2t_chr${CHR}_features_summary.json"

        if [[ ! -f "${G_VCF}" ]]; then
            log "WARNING: Missing GRCh38 chr${CHR} VCF — skipping chromosome. (${G_VCF})"
            continue
        fi
        if [[ ! -f "${T_VCF}" ]]; then
            log "WARNING: Missing T2T chr${CHR} VCF — skipping chromosome. (${T_VCF})"
            continue
        fi

        log "  Extracting chr${CHR} features (GRCh38)..."
        "${PY}" "${CODE_DIR}/02b_extract_variant_features.py" \
            --vcf "${G_VCF}" \
            --metadata "${META_CLEAN}" \
            --reference-label GRCh38 \
            --max-variants "${VARIANTS_PER_CHR}" \
            --maf-threshold "${MAF}" \
            --call-rate-threshold "${CALL_RATE}" \
            --out-features-csv "${G_OUT}" \
            --out-summary-json "${G_SUM}" \
            >> "${RUN_LOG}" 2>&1

        log "  Extracting chr${CHR} features (T2T-CHM13)..."
        "${PY}" "${CODE_DIR}/02b_extract_variant_features.py" \
            --vcf "${T_VCF}" \
            --metadata "${META_CLEAN}" \
            --reference-label T2T_CHM13 \
            --max-variants "${VARIANTS_PER_CHR}" \
            --maf-threshold "${MAF}" \
            --call-rate-threshold "${CALL_RATE}" \
            --out-features-csv "${T_OUT}" \
            --out-summary-json "${T_SUM}" \
            >> "${RUN_LOG}" 2>&1

        G_CHR_CSVS+=("${G_OUT}")
        T_CHR_CSVS+=("${T_OUT}")
    done

    if [[ ${#G_CHR_CSVS[@]} -eq 0 ]]; then
        log "ERROR: No chromosomes extracted successfully. Aborting."
        exit 1
    fi

    log "  Merging ${#G_CHR_CSVS[@]} chromosome CSVs into genome-wide GRCh38 matrix..."
    "${PY}" "${CODE_DIR}/02f_merge_chromosome_features.py" \
        --input-csvs "${G_CHR_CSVS[@]}" \
        --out-features-csv "${GRCH38_MERGED}" \
        --out-summary-json "${MANIFESTS_DIR}/grch38_merge_summary.json" \
        | tee -a "${RUN_LOG}"

    log "  Merging ${#T_CHR_CSVS[@]} chromosome CSVs into genome-wide T2T matrix..."
    "${PY}" "${CODE_DIR}/02f_merge_chromosome_features.py" \
        --input-csvs "${T_CHR_CSVS[@]}" \
        --out-features-csv "${T2T_MERGED}" \
        --out-summary-json "${MANIFESTS_DIR}/t2t_merge_summary.json" \
        | tee -a "${RUN_LOG}"

    log "  Aligning common variants between references..."
    "${PY}" "${CODE_DIR}/02c_align_reference_features.py" \
        --grch38-features "${GRCH38_MERGED}" \
        --t2t-features "${T2T_MERGED}" \
        --min-common-patients 100 \
        --min-features-per-ref 20 \
        --out-grch38-aligned "${GRCH38_ALIGNED}" \
        --out-t2t-aligned "${T2T_ALIGNED}" \
        --out-summary-json "${ALIGN_SUMMARY}" \
        | tee -a "${RUN_LOG}"

    log "Phase 2 complete."
else
    log "--- Phase 2 skipped (--skip-phase2) ---"
    for f in "${GRCH38_ALIGNED}" "${T2T_ALIGNED}"; do
        if [[ ! -f "${f}" ]]; then
            log "ERROR: --skip-phase2 set but required file is missing: ${f}"; exit 1
        fi
    done
fi

# -----------------------------------------------------------------------
# PHASE 3: Training and evaluation — all models, both references
# -----------------------------------------------------------------------
log "--- Phase 3: Model training and sex-stratified evaluation ---"

ANALYSIS_DIR="${METRICS_DIR}/genome_wide"
PRED_ANALYSIS_DIR="${PREDS_DIR}/genome_wide"
mkdir -p "${ANALYSIS_DIR}" "${PRED_ANALYSIS_DIR}"

for MODEL in ${MODELS}; do
    log "  Training: ${MODEL} on GRCh38..."
    "${PY}" "${CODE_DIR}/03_train_models.py" \
        --features-csv "${GRCH38_ALIGNED}" \
        --split-manifest "${SPLIT_CSV}" \
        --model "${MODEL}" \
        --seed "${SEED}" \
        --max-repeats "${REPEATS}" \
        --metrics-out "${ANALYSIS_DIR}/grch38_${MODEL}_metrics.csv" \
        --predictions-out "${PRED_ANALYSIS_DIR}/grch38_${MODEL}_predictions.csv" \
        --summary-out "${ANALYSIS_DIR}/grch38_${MODEL}_summary.json" \
        | tee -a "${RUN_LOG}"

    log "  Training: ${MODEL} on T2T-CHM13..."
    "${PY}" "${CODE_DIR}/03_train_models.py" \
        --features-csv "${T2T_ALIGNED}" \
        --split-manifest "${SPLIT_CSV}" \
        --model "${MODEL}" \
        --seed "${SEED}" \
        --max-repeats "${REPEATS}" \
        --metrics-out "${ANALYSIS_DIR}/t2t_${MODEL}_metrics.csv" \
        --predictions-out "${PRED_ANALYSIS_DIR}/t2t_${MODEL}_predictions.csv" \
        --summary-out "${ANALYSIS_DIR}/t2t_${MODEL}_summary.json" \
        | tee -a "${RUN_LOG}"
done

log "Phase 3 complete."

# -----------------------------------------------------------------------
# PHASE 4: Primary estimand analysis (ΔΔ) for each model
# -----------------------------------------------------------------------
log "--- Phase 4: ΔΔ analysis ---"

FINAL_DIR="${METRICS_DIR}/final"
mkdir -p "${FINAL_DIR}"

for MODEL in ${MODELS}; do
    G_METRICS="${ANALYSIS_DIR}/grch38_${MODEL}_metrics.csv"
    T_METRICS="${ANALYSIS_DIR}/t2t_${MODEL}_metrics.csv"

    if [[ ! -f "${G_METRICS}" || ! -f "${T_METRICS}" ]]; then
        log "  WARNING: Metrics missing for ${MODEL} — skipping ΔΔ analysis."
        continue
    fi

    log "  Analyzing ΔΔ for model: ${MODEL}..."
    "${PY}" "${CODE_DIR}/04_analyze_results.py" \
        --grch38-metrics "${G_METRICS}" \
        --t2t-metrics "${T_METRICS}" \
        --bootstrap-iters 5000 \
        --permutation-iters 5000 \
        --seed "${SEED}" \
        --out-csv "${FINAL_DIR}/${MODEL}_delta_delta_by_repeat.csv" \
        --out-json "${FINAL_DIR}/${MODEL}_delta_delta_summary.json" \
        | tee -a "${RUN_LOG}"
done

log "Phase 4 complete."

# -----------------------------------------------------------------------
# PHASE 5: Data manifest / reproducibility snapshot
# -----------------------------------------------------------------------
log "--- Phase 5: Generating reproducibility manifest ---"
bash "${CODE_DIR}/00_generate_data_manifest.sh" 2>&1 | tee -a "${RUN_LOG}"
log "Phase 5 complete."

log ""
log "=== Pipeline complete ==="
log "Results: ${RESULTS_DIR}"
log "Run log: ${RUN_LOG}"
log "Primary outputs:"
for MODEL in ${MODELS}; do
    SUMMARY="${FINAL_DIR}/${MODEL}_delta_delta_summary.json"
    if [[ -f "${SUMMARY}" ]]; then
        log "  [${MODEL}] $(cat ${SUMMARY} | "${PY}" -c "import sys,json; d=json.load(sys.stdin); print(f'ΔΔ={d[\"delta_delta_mean\"]:.4f}, CI=[{d[\"delta_delta_bootstrap_ci_95\"][0]:.4f},{d[\"delta_delta_bootstrap_ci_95\"][1]:.4f}], p={d[\"paired_sign_permutation_p_value\"]:.4f}, supports_H1={d[\"direction_supports_hypothesis\"]}')")"
    fi
done
