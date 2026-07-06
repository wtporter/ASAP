#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
R_ENV="${REPO_ROOT}/nextflow/work/r_test_env"
RSCRIPT="${R_ENV}/bin/Rscript"
TEST_DIR="${REPO_ROOT}/asap_tools/tests/testthat"
LOG_DIR="${REPO_ROOT}/asap_tools/tests/logs"
mkdir -p "${LOG_DIR}"

if [[ ! -x "${RSCRIPT}" ]]; then
  echo "R test env not found at ${R_ENV} — building from r_env.yml..."
  mamba env create \
    -f "${REPO_ROOT}/nextflow/modules/asap_tools/r_env.yml" \
    -p "${R_ENV}"
  echo "Env built successfully."
fi

sbatch \
  --partition=compute \
  --cpus-per-task=4 \
  --mem=16G \
  --time=0:30:00 \
  --job-name=asap_r_tests \
  --output="${LOG_DIR}/asap_r_tests_%j.log" \
  --wrap="${RSCRIPT} -e \"testthat::test_dir('${TEST_DIR}', reporter='progress')\""

echo "Submitted. Logs will appear as ${LOG_DIR}/asap_r_tests_<JOBID>.log"
