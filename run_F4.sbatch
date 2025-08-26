#!/bin/bash
#SBATCH --job-name=gb_F4
#SBATCH --partition=turing-long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --chdir=/home/fspreafichi/julia_grobner_pipeline
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

mkdir -p logs results

if [[ $# -lt 1 ]]; then
  echo "Usage: sbatch run_F4.sbatch <input_file.in> [threads]"
  exit 1
fi

INPUT_FILE="$1"
THREADS="${2:-${SLURM_CPUS_PER_TASK:-8}}"

if [[ ! -f "$INPUT_FILE" ]]; then
  echo "Input file not found: $INPUT_FILE"
  exit 2
fi

# Avoid BLAS/OpenMP storms; set Julia threads
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export JULIA_NUM_THREADS="$THREADS"

echo "[$(date)] Node(s): ${SLURM_JOB_NODELIST:-unknown}"
echo "[$(date)] CWD: $(pwd)"
echo "[$(date)] Input: $INPUT_FILE"
echo "[$(date)] Threads: $THREADS"

module load julia
which julia
julia -v
julia -e 'println("Julia smoke test ok; Threads=", Threads.nthreads())'

julia scripts/solve_F4_from_file.jl "$INPUT_FILE" "$THREADS"

echo "[$(date)] Listing results/:"
ls -lh results || true
echo "[$(date)] ✅ Done"

