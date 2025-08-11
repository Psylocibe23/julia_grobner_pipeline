#!/bin/bash
#SBATCH --job-name=fglm
#SBATCH --output=logs/%x_%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem=32G

# Usage:
#   sbatch run_FGLM.sbatch <F4_result_file.txt>

set -euo pipefail

INPUT_FILE="${1:-}"
if [[ -z "${INPUT_FILE}" ]]; then
  echo "Usage: sbatch run_FGLM.sbatch <F4_result_file.txt>"
  exit 1
fi
if [[ ! -f "${INPUT_FILE}" ]]; then
  echo "ERROR: file not found: ${INPUT_FILE}"
  exit 1
fi

module load sagemath/10.6-singularity
sage scripts/convert_to_lex_fglm.sage "${INPUT_FILE}"
