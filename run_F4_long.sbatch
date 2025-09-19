#!/bin/bash
#SBATCH --job-name=gb_F4_long
#SBATCH --qos=tesi_spreafichi
#SBATCH --partition=turing-wide
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=960G                 
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail
mkdir -p logs results

if [[ $# -lt 1 ]]; then
  echo "Usage: sbatch run_F4_long.sbatch <input_file.in> [threads]"
  exit 1
fi

INPUT_FILE="$1"
THREADS="${2:-${SLURM_CPUS_PER_TASK:-32}}"

[[ -f "$INPUT_FILE" ]] || { echo "Input not found: $INPUT_FILE"; exit 2; }

# avoid oversubscription from math libs (F4 uses its own threads)
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export JULIA_NUM_THREADS=1

echo "[$(date)] Node(s): ${SLURM_JOB_NODELIST:-unknown}"
echo "[$(date)] Input: $INPUT_FILE"
echo "[$(date)] Threads to F4 (nr_thrds): $THREADS"

# modules
if ! type module >/dev/null 2>&1; then
  [[ -f /etc/profile ]] && source /etc/profile
  [[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh
fi
module load julia

which julia
julia -v

# run
julia scripts/solve_F4_from_file_parsing.jl "$INPUT_FILE" "$THREADS"

echo "[$(date)] Listing results/:"
ls -lh results || true
echo "[$(date)] ✅ Done"
