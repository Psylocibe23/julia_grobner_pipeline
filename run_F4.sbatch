#!/bin/bash -l
#SBATCH --job-name=gb_F4
#SBATCH --partition=turing-long
#SBATCH --chdir=/home/fspreafichi/julia_grobner_pipeline
#SBATCH --output=gb_F4_%j.out         # write to CWD (guaranteed to exist)
#SBATCH --error=gb_F4_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00

# Usage:
#   sbatch run_F4.sbatch data/HFE_n5_D6.in [threads]

set -euo pipefail
set -x

# Make places for *your* logs/results (this happens after Slurm has opened stdout/err)
mkdir -p logs results

if [[ $# -lt 1 ]]; then
  echo "Usage: sbatch run_F4.sbatch <input_file.in> [threads]"
  exit 1
fi

INPUT_FILE="$1"
if [[ ! -f "$INPUT_FILE" ]]; then
  echo "ERROR: file not found: $INPUT_FILE"
  exit 2
fi

THREADS="${2:-${SLURM_CPUS_PER_TASK:-8}}"

# keep math libs quiet; give Julia threads
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export JULIA_NUM_THREADS="$THREADS"

echo "[$(date)] Node(s): ${SLURM_JOB_NODELIST:-unknown}"
echo "[$(date)] CWD: $(pwd)"
echo "[$(date)] Input: $INPUT_FILE"
echo "[$(date)] Threads: $THREADS"

module load julia || { echo "module load julia failed"; exit 3; }
which julia
julia -v
julia -e 'println("Julia smoke test ok; Threads=",Threads.nthreads())'

julia scripts/solve_F4_from_file.jl "$INPUT_FILE" "$THREADS"

echo "[$(date)] Listing results/:"
ls -lh results || true
