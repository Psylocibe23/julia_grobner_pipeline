#!/bin/bash
#SBATCH --job-name=genhfe
#SBATCH --output=logs/genhfe_%j.out
#SBATCH --error=logs/genhfe_%j.err
#SBATCH --partition=turing-long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32          # the script will auto-limit workers as needed
#SBATCH --mem=128G
#SBATCH --time=08:00:00

set -euo pipefail

module load sagemath/10.6-singularity

# -----------------------
# CLI parameters (all optional after N and D)
# Usage:
#   sbatch run_generate_HFE_parallel.sbatch N D [SEED] [PROB_QUAD] [PROB_LIN] [RANK_MIN] [ALLOW_FB] [MAX_MAPS] [MAX_TRIES] [WORKERS] [FORCE]
# Examples:
#   sbatch run_generate_HFE_parallel.sbatch 80 96
#   sbatch run_generate_HFE_parallel.sbatch 80 96 123 0.65 0.60 - true 30 4096 auto true
# -----------------------
N=${1:-80}
D=${2:-96}
SEED=${3:-0}
PROB_QUAD=${4:-0.65}
PROB_LIN=${5:-0.60}
RANK_MIN=${6:--}       # "-" = auto (n-1)
ALLOW_FB=${7:-true}
MAX_MAPS=${8:-30}
MAX_TRIES=${9:-4096}
WORKERS=${10:-auto}    # "auto" → let the script decide; integer to override
FORCE=${11:-false}

mkdir -p logs data/hfe_instances

# Avoid BLAS/OMP oversubscription: our script uses process parallelism.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# Force overwrite via env (preferred so we don't disturb positional args)
if [[ "$FORCE" == "true" ]]; then
  export HFE_FORCE=1
fi

# Let the Sage script choose the default outfile path
OUTFILE="-"

# If WORKERS is "auto", omit it so the script auto-picks from SLURM_CPUS_PER_TASK;
# otherwise pass the integer you gave. We don't pass a final 'force' arg—env above handles it.
if [[ "$WORKERS" == "auto" ]]; then
  echo "[$(date)] Node(s): $SLURM_JOB_NODELIST"
  echo "[$(date)] Generating HFE($N,$D) with up to $SLURM_CPUS_PER_TASK CPUs; workers=auto"
  sage scripts/construct_HFE_parallelized.sage \
      "$N" "$D" "$OUTFILE" "$SEED" "$PROB_QUAD" "$PROB_LIN" "$RANK_MIN" "$ALLOW_FB" \
      "$MAX_MAPS" "$MAX_TRIES"
else
  echo "[$(date)] Node(s): $SLURM_JOB_NODELIST"
  echo "[$(date)] Generating HFE($N,$D) with up to $SLURM_CPUS_PER_TASK CPUs; workers=$WORKERS"
  sage scripts/construct_HFE_parallelized.sage \
      "$N" "$D" "$OUTFILE" "$SEED" "$PROB_QUAD" "$PROB_LIN" "$RANK_MIN" "$ALLOW_FB" \
      "$MAX_MAPS" "$MAX_TRIES" "$WORKERS"
fi

echo "[$(date)] ✅ Done"
