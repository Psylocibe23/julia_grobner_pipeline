#!/bin/bash
#SBATCH --job-name=gb_F4
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --chdir=/home/fspreafichi/julia_grobner_pipeline
#SBATCH --ntasks=1
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
set -x

if [[ $# -lt 1 ]]; then
  echo "Usage: sbatch run_F4.sbatch <input_file.in> [threads]"
  exit 1
fi

INPUT_FILE="$1"
if [[ ! -f "$INPUT_FILE" ]]; then
  echo "ERROR: file not found: $INPUT_FILE"
  exit 1
fi

# Decide threads: arg > SLURM_CPUS_PER_TASK > 8
THREADS="${2:-${SLURM_CPUS_PER_TASK:-8}}"

# Make sure relative output paths exist
mkdir -p logs results

# Avoid BLAS/OpenMP oversubscription; enable Julia threads
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export JULIA_NUM_THREADS="$THREADS"

echo "[$(date)] Node(s): ${SLURM_JOB_NODELIST:-unknown}"
echo "[$(date)] CWD: $(pwd)"
echo "[$(date)] Input: $INPUT_FILE"
echo "[$(date)] Threads: $THREADS"

module load julia
julia scripts/solve_F4_from_file.jl "$INPUT_FILE" "$THREADS"

echo "[$(date)] Listing results/ after run:"
ls -lh results || true
