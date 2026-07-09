#!/usr/bin/env bash
#SBATCH --job-name=ASAP_prebuild_envs
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=logs/nf-test/prebuild_envs_%j.out
#
# Prebuild the two conda environments the ASAP pipeline uses, ONCE, into a
# shared cache directory. run_tests.sh then points every (parallel) test job at
# these prebuilt envs via ASAP_CONDA_ENV / R_CONDA_ENV, so no job has to build
# an env itself. This eliminates the mamba lock contention on the shared pkgs
# cache ("Could not set lock ... /tgen_labs/EPIC/miniconda3/pkgs") and the
# concurrent-conda corruption races ("NoSuchFileException ... *.pyc") that fail
# tests when many jobs build the same env at the same time.
#
# Only two envs are conda-based; all other modules (FastQC, MultiQC, iVAR,
# aligners) use Singularity containers and are unaffected.
#
# Usage:
#   sbatch prebuild_envs.sh                  # build (or reuse) both envs
#   ./prebuild_envs.sh                       # run directly (login/interactive node)
#   ASAP_ENV_CACHE=/path ./prebuild_envs.sh  # override the cache location
#   ASAP_FORCE_REBUILD=1 ./prebuild_envs.sh  # remove and rebuild existing envs
#
# The cache defaults to nextflow/work/conda_envs (gitignored). Building is
# idempotent: an env that already exists is left in place unless
# ASAP_FORCE_REBUILD=1.

set -euo pipefail

source /tgen_labs/EPIC/miniconda3/etc/profile.d/conda.sh
conda activate ASAP_nextflow_env

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "${SCRIPT_DIR}"
mkdir -p logs/nf-test

ASAP_ENV_CACHE="${ASAP_ENV_CACHE:-${SCRIPT_DIR}/work/conda_envs}"
mkdir -p "${ASAP_ENV_CACHE}"

# env name -> source YAML
declare -A ENVS=(
    ["asap_env"]="modules/asap/asap_env.yml"
    ["r_env"]="modules/asap_tools/r_env.yml"
)

build_env() {
    local name="$1" yml="$2"
    local prefix="${ASAP_ENV_CACHE}/${name}"

    if [[ "${ASAP_FORCE_REBUILD:-0}" == "1" && -d "${prefix}" ]]; then
        echo ">> Removing existing env for rebuild: ${prefix}"
        rm -rf "${prefix}"
    fi

    if [[ -d "${prefix}" ]]; then
        echo ">> Reusing existing env: ${prefix}"
        return 0
    fi

    echo ">> Building ${name} from ${yml} -> ${prefix}"
    # Build directly into the final prefix. Conda envs are NOT relocatable
    # (shebangs and activation scripts hardcode the build prefix), so we must
    # not build elsewhere and move. If the build fails or is interrupted, remove
    # the partial prefix so a later run doesn't mistake it for a complete env.
    trap 'echo ">> Build failed; removing partial env ${prefix}"; rm -rf "${prefix}"' ERR
    mamba env create --yes --prefix "${prefix}" --file "${yml}"
    trap - ERR
    echo ">> Done: ${prefix}"
}

echo "ASAP conda env cache: ${ASAP_ENV_CACHE}"
for name in "${!ENVS[@]}"; do
    build_env "${name}" "${ENVS[$name]}"
done

echo ""
echo "All envs ready under ${ASAP_ENV_CACHE}:"
for name in "${!ENVS[@]}"; do
    echo "  ${name} -> ${ASAP_ENV_CACHE}/${name}"
done
