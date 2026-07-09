#!/usr/bin/env bash
set -euo pipefail

# Runs the asap_tools R unit/integration tests (testthat).
#
# Env resolution (first match wins):
#   1. $ASAP_R_ENV                     — explicit conda env prefix, if valid
#   2. any env under $NF_CONDA_DIR     — prebuilt Nextflow conda envs
#      (default /scratch/tporter/nf_envs/conda) and the repo's
#      nextflow/work/conda that has Rscript + genbankr + testthat
#   3. mamba create from r_env.yml     — last resort, into $ASAP_R_ENV
#
# A valid env must have bin/Rscript AND lib/R/etc/ldpaths (a missing ldpaths is
# the failure mode of a half-built env — its Rscript aborts on startup).
#
# Run mode: local by default; pass --sbatch to submit to SLURM.

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
TEST_DIR="${REPO_ROOT}/asap_tools/tests/testthat"
LOG_DIR="${REPO_ROOT}/asap_tools/tests/logs"
NF_CONDA_DIR="${NF_CONDA_DIR:-/scratch/tporter/nf_envs/conda}"
R_ENV_FALLBACK="${ASAP_R_ENV:-${NF_CONDA_DIR}/asap_r_test_env}"
mkdir -p "${LOG_DIR}"

SUBMIT_SBATCH=0
[[ "${1:-}" == "--sbatch" ]] && SUBMIT_SBATCH=1

# True if $1 is a usable R env prefix (has Rscript and ldpaths).
is_valid_env() {
  [[ -x "$1/bin/Rscript" && -f "$1/lib/R/etc/ldpaths" ]]
}

# True if $1's R library has genbankr and testthat.
has_r_deps() {
  [[ -d "$1/lib/R/library/genbankr" && -d "$1/lib/R/library/testthat" ]]
}

RSCRIPT=""

# 1. explicit override
if [[ -n "${ASAP_R_ENV:-}" ]] && is_valid_env "${ASAP_R_ENV}"; then
  RSCRIPT="${ASAP_R_ENV}/bin/Rscript"
fi

# 2. discover a prebuilt env
if [[ -z "${RSCRIPT}" ]]; then
  for base in "${NF_CONDA_DIR}" "${REPO_ROOT}/nextflow/work/conda"; do
    [[ -d "${base}" ]] || continue
    for env in "${base}"/*/; do
      env="${env%/}"
      if is_valid_env "${env}" && has_r_deps "${env}"; then
        RSCRIPT="${env}/bin/Rscript"
        break 2
      fi
    done
  done
fi

# 3. build from spec as a last resort
if [[ -z "${RSCRIPT}" ]]; then
  echo "No prebuilt R env found; creating ${R_ENV_FALLBACK} from r_env.yml..."
  mamba env create --yes \
      -f "${REPO_ROOT}/nextflow/modules/asap_tools/r_env.yml" \
      -p "${R_ENV_FALLBACK}" 2>/dev/null || \
  mamba env update --prune --yes \
      -f "${REPO_ROOT}/nextflow/modules/asap_tools/r_env.yml" \
      -p "${R_ENV_FALLBACK}"
  is_valid_env "${R_ENV_FALLBACK}" || { echo "ERROR: built env is unusable (no ldpaths)"; exit 1; }
  RSCRIPT="${R_ENV_FALLBACK}/bin/Rscript"
fi

echo "Using Rscript: ${RSCRIPT}"

R_CMD="testthat::test_dir('${TEST_DIR}', reporter='summary', stop_on_failure=TRUE)"

if [[ "${SUBMIT_SBATCH}" -eq 1 ]]; then
  sbatch \
    --partition=compute \
    --cpus-per-task=4 \
    --mem=16G \
    --time=0:30:00 \
    --job-name=asap_r_tests \
    --output="${LOG_DIR}/asap_r_tests_%j.log" \
    --wrap="${RSCRIPT} -e \"${R_CMD}\""
  echo "Submitted. Logs will appear as ${LOG_DIR}/asap_r_tests_<JOBID>.log"
else
  "${RSCRIPT}" -e "${R_CMD}"
fi
