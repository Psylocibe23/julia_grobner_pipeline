#!/bin/bash
#SBATCH --job-name=sysdiag
#SBATCH --partition=turing-long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=03:00:00
# Write logs next to the job so no missing-dir issues:
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
# Work inside your repo so relative paths like data/... resolve:
#SBATCH --chdir=/home/fspreafichi/julia_grobner_pipeline

set -euo pipefail
set -x  # echo commands for debugging

module load sagemath/10.6-singularity

# Keep math libs from oversubscribing threads
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# The Sage script writes its own log to logs/... so make sure it exists
mkdir -p logs

if [[ $# -lt 1 ]]; then
  echo "Usage: sbatch run_system_diagnosis.sbatch <input.in> [--dim=SECONDS] [--hilbert=SECONDS]"
  exit 1
fi

INPUT_FILE="$1"
shift || true

if [[ ! -f "$INPUT_FILE" ]]; then
  echo "Input file not found: $INPUT_FILE"
  exit 2
fi

echo "[$(date)] Node(s): ${SLURM_JOB_NODELIST:-unknown}"
echo "[$(date)] Running diagnosis for: $INPUT_FILE"
echo "[$(date)] Extra args → $*"

# Run it. Any flags after the file (e.g. --dim=1800 --hilbert=600) are forwarded.
sage scripts/system_diagnosis.sage "$INPUT_FILE" "$@"

echo "[$(date)] ✅ Done"
