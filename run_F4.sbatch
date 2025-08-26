#!/bin/bash
#SBATCH --job-name=gb_F4
# (Make sure the logs/ directory exists before submitting, or use the plain form below)
#SBATCH --output=logs/%x_%j.out
#SBATCH --ntasks=1
# Slurm Directives MUST be literal constants; override at submit time if needed.
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --mem=32G

# Usage:
#   sbatch run_F4.sbatch <input_file.in> [threads]
# Examples:
#   sbatch run_F4.sbatch data/HFE_n5_D6.in
#   sbatch --cpus-per-task=16 run_F4.sbatch data/HFE_n32_D64.in
#   sbatch run_F4.sbatch data/HFE_n32_D64.in 12

set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: sbatch run_F4.sbatch <input_file.in> [threads]"
  exit 1
fi

INPUT_FILE="$1"
if [[ ! -f "$INPUT_FILE" ]]; then
  echo "ERROR: file not found: $INPUT_FILE"
  exit 1
fi

# Decide how many threads to use:
# 1) if a positional [threads] was given, use it,
# 2) else, use SLURM_CPUS_PER_TASK (e.g. when you submit with --cpus-per-task=...),
# 3) else, fall back to 8 (matches the #SBATCH default).
THREADS_DEFAULT="${SLURM_CPUS_PER_TASK:-8}"
THREADS="${2:-$THREADS_DEFAULT}"

# Keep math libs from oversubscribing threads (Singular/BLAS)
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
# If your Julia script uses Julia threads, this enables them:
export JULIA_NUM_THREADS="$THREADS"

echo "[$(date)] Node(s): ${SLURM_JOB_NODELIST:-unknown}"
echo "[$(date)] Using THREADS=$THREADS  (SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-unset})"
echo "[$(date)] Input: $INPUT_FILE"

module load julia
julia scripts/solve_F4_from_file.jl "$INPUT_FILE" "$THREADS"
