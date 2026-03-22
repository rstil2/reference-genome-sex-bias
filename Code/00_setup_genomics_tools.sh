#!/usr/bin/env bash
set -euo pipefail

CONDA_BIN="${CONDA_BIN:-/Users/stillwell/miniconda3/bin/conda}"
CONDA_PREFIX_TARGET="${CONDA_PREFIX_TARGET:-/Users/stillwell/miniconda3}"

echo "Installing genomics CLI dependencies into: ${CONDA_PREFIX_TARGET}"
set +e
"${CONDA_BIN}" install -y -p "${CONDA_PREFIX_TARGET}" -c conda-forge -c bioconda \
  bcftools htslib plink2
status=$?
set -e

if [[ ${status} -ne 0 ]]; then
  echo "plink2 not available in current channels; falling back to plink."
  "${CONDA_BIN}" install -y -p "${CONDA_PREFIX_TARGET}" -c conda-forge -c bioconda \
    bcftools htslib plink
fi

echo "Dependency install complete."
echo "Versions:"
"${CONDA_PREFIX_TARGET}/bin/bcftools" --version | head -n 1
if [[ -x "${CONDA_PREFIX_TARGET}/bin/plink2" ]]; then
  "${CONDA_PREFIX_TARGET}/bin/plink2" --version | head -n 1
elif [[ -x "${CONDA_PREFIX_TARGET}/bin/plink" ]]; then
  "${CONDA_PREFIX_TARGET}/bin/plink" --version | head -n 1
else
  echo "WARNING: Neither plink2 nor plink found on expected path."
fi
