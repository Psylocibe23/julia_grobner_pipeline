#!/bin/bash
#SBATCH --job-name=gb_F4_long
#SBATCH --qos=tesi_spreafichi
#SBATCH --partition=turing-wide,Main,lovelace
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=512G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail
mkdir -p logs results

if [[ $# -lt 1 ]]; then
  echo "Usage: sbatch run_F4_long.sbatch <input_file.in> [threads]"
  exit 1
fi

INPUT_FILE="$1"
THREADS="${2:-${SLURM_CPUS_PER_TASK:-8}}"

[[ -f "$INPUT_FILE" ]] || { echo "Input not found: $INPUT_FILE"; exit 2; }

# Avoid oversubscription from math libs
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export JULIA_NUM_THREADS="$THREADS"

echo "[$(date)] Node(s): ${SLURM_JOB_NODELIST:-unknown}"
echo "[$(date)] CWD: $(pwd)"
echo "[$(date)] Input: $INPUT_FILE"
echo "[$(date)] Threads: $THREADS"

# Load Julia
if ! type module >/dev/null 2>&1; then
  [[ -f /etc/profile ]] && source /etc/profile
  [[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh
fi
module load julia

which julia
julia -v
julia -e 'println("Julia ok; Threads=", Threads.nthreads())'

# Optional: use local depot/project to reduce network pressure (uncomment if helpful)
# export JULIA_DEPOT_PATH="$PWD/.julia_depot"
# export JULIA_PROJECT="$PWD"
# julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

# Run
julia scripts/solve_F4_from_file_parsing.jl "$INPUT_FILE" "$THREADS"

echo "[$(date)] Listing results/:"
ls -lh results || true
echo "[$(date)] ✅ Done"
