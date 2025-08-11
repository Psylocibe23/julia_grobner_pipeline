#!/bin/bash
#SBATCH --job-name=gb_F4
#SBATCH --output=logs/%x_%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${2:-8}
#SBATCH --time=01:00:00
#SBATCH --mem=32G

# Usage:
#   sbatch run_F4.sbatch <input_file.in> [threads]

set -euo pipefail

INPUT_FILE="${1:-}"
THREADS="${2:-8}"

if [[ -z "${INPUT_FILE}" ]]; then
  echo "Usage: sbatch run_F4.sbatch <input_file.in> [threads]"
  exit 1
fi
if [[ ! -f "${INPUT_FILE}" ]]; then
  echo "ERROR: file not found: ${INPUT_FILE}"
  exit 1
fi

module load julia
julia scripts/solve_F4_from_file.jl "${INPUT_FILE}" "${THREADS}"
