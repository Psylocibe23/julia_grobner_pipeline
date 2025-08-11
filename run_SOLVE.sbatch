#!/bin/bash
#SBATCH --job-name=solve
#SBATCH --output=logs/%x_%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem=32G

# Usage:
#   sbatch run_SOLVE.sbatch <LEX_basis_file.txt>

set -euo pipefail

INPUT_FILE="${1:-}"
if [[ -z "${INPUT_FILE}" ]]; then
  echo "Usage: sbatch run_SOLVE.sbatch <LEX_basis_file.txt>"
  exit 1
fi
if [[ ! -f "${INPUT_FILE}" ]]; then
  echo "ERROR: file not found: ${INPUT_FILE}"
  exit 1
fi

module load sagemath/10.6-singularity
sage scripts/extract_solutions_from_lex.sage "${INPUT_FILE}"
