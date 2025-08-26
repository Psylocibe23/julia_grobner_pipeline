#!/usr/bin/env bash
#SBATCH --job-name=gb_F4
#SBATCH --partition=turing-long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --chdir=/home/fspreafichi/julia_grobner_pipeline
#SBATCH --output=%x_%j.out        # write logs in CWD so no missing-dir issue
#SBATCH --error=%x_%j.err

set -euo pipefail
set -x

# Args: sbatch run_F4.sbatch data/HFE_n6_D6.in [threads]
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

# Keep math libs from oversubscribing; give Julia threads
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export JULIA_NUM_THREADS="$THREADS"

echo "[$(date)] Node(s): ${SLURM_JOB_NODELIST:-unknown}"
echo "[$(date)] CWD: $(pwd)"
echo "[$(date)] Input: $INPUT_FILE"
echo "[$(date)] Threads: $THREADS"

# Make sure module command exists even on non-login shells
if ! type module >/dev/null 2>&1; then
  # Common init locations; harmless if missing
  [[ -f /etc/profile ]] && source /etc/profile
  [[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh
fi

module load julia
which julia
julia -v
julia -e 'println("Julia smoke test ok; Threads=", Threads.nthreads())'

# Ensure your output folders exist (for files created by the Julia script)
mkdir -p logs results

# Run the solver
julia scripts/solve_F4_from_file.jl "$INPUT_FILE" "$THREADS"

echo "[$(date)] Listing results/:"
ls -lh results || true
echo "[$(date)] ✅ Done"
