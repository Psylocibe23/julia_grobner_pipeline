#!/bin/bash
#SBATCH --job-name=gb_F5
#SBATCH --partition=turing-long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1          # F5 path is effectively single-threaded
#SBATCH --mem=64G                  
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

mkdir -p logs results

# Usage/args (same style as F4)
if [[ $# -lt 1 ]]; then
  echo "Usage: sbatch run_F5.sbatch <input_file.in>"
  exit 1
fi

INPUT_FILE="$1"
if [[ ! -f "$INPUT_FILE" ]]; then
  echo "Input file not found: $INPUT_FILE"
  exit 2
fi

# Pin all threading to 1: F5 itself is single-threaded; keep BLAS quiet too
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export JULIA_NUM_THREADS=1

echo "[$(date)] Node(s): ${SLURM_JOB_NODELIST:-unknown}"
echo "[$(date)] CWD: $(pwd)"
echo "[$(date)] Input: $INPUT_FILE"
echo "[$(date)] Threads pinned: $JULIA_NUM_THREADS (F5 single-threaded wrapper)"

# Make 'module' available in non-login shells, then load Julia
if ! type module >/dev/null 2>&1; then
  [[ -f /etc/profile ]] && source /etc/profile
  [[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh
fi
module load julia

which julia
julia -v
julia -e 'println("Julia smoke test ok; Threads=", Threads.nthreads())'

# Run the F5 solver (writes logs/<...>.log and results/<...>.txt)
julia scripts/solve_F5_from_file.jl "$INPUT_FILE" 

echo "[$(date)] Listing results/:"
ls -lh results || true
echo "[$(date)] ✅ Done"
