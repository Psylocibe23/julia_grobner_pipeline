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

# Match sysdiag: make sure these exist where Slurm writes logs
mkdir -p logs results

# Args: identical style to sysdiag
if [[ $# -lt 1 ]]; then
  echo "Usage: sbatch run_F4.sbatch <input_file.in> [threads]"
  exit 1
fi

INPUT_FILE="$1"
shift || true
if [[ ! -f "$INPUT_FILE" ]]; then
  echo "Input file not found: $INPUT_FILE"
  exit 2
fi

# Threads: [arg] > SLURM_CPUS_PER_TASK > 8
THREADS="${1:-${SLURM_CPUS_PER_TASK:-8}}"

# Keep math libs from oversubscribing
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export JULIA_NUM_THREADS="$THREADS"

echo "[$(date)] Node(s): ${SLURM_JOB_NODELIST:-unknown}"
echo "[$(date)] CWD: $(pwd)"
echo "[$(date)] Input: $INPUT_FILE"
echo "[$(date)] Threads: $THREADS"

# Clean module env (optional but helps) then load Julia
module purge >/dev/null 2>&1 || true
module load julia

# Smoke test so we see immediate failures in this job’s stdout
which julia || { echo "julia not in PATH"; exit 3; }
julia -v
julia -e 'println("Julia smoke test ok; Threads=", Threads.nthreads())'

# If AlgebraicSolving isn’t preinstalled, this will error *and* the message
# will land in logs/gb_F4_<jobid>.err/out so you can see it.
julia scripts/solve_F4_from_file.jl "$INPUT_FILE" "$THREADS"

echo "[$(date)] Listing results/:"
ls -lh results || true
echo "[$(date)] ✅ Done"
