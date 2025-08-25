#!/bin/bash
#SBATCH --job-name=sysdiag
#SBATCH --partition=turing-long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

# Ensure logs dir exists (matches --output/--error paths)
mkdir -p logs

# Keep math libs from oversubscribing threads
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

module load sagemath/10.6-singularity

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
if [[ $# -gt 0 ]]; then
  echo "[$(date)] Extra args → $*"
fi

# Run the Sage script; any flags are forwarded (e.g. --dim=1800 --hilbert=600)
sage scripts/system_diagnosis.sage "$INPUT_FILE" "$@"

echo "[$(date)] ✅ Done"
