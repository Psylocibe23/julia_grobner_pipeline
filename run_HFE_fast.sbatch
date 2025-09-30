#!/bin/bash
#SBATCH --job-name=hfe_fast
#SBATCH --output=logs/hfe_fast_%j.out
#SBATCH --error=logs/hfe_fast_%j.err
#SBATCH --qos=tesi_spreafichi
#SBATCH --partition=turing-wide
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1        # single-threaded generator; keep tiny to avoid QOS limits
#SBATCH --mem=16G                # plenty for generation
#SBATCH --time=02:00:00          # walltime > internal 90-min budget for clean wrap-up

set -euo pipefail

# be explicit: no hidden multithreading
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

module load sagemath/10.6-singularity

# Usage:
#   sbatch run_HFE_fast.sbatch N D [SEED|-] [TIME_BUDGET_SEC] [MAX_MAPS] [SECRET_TRIES] [RANK_MIN|-] [FORCE(true/false)]
#
# Behavior:
#   - The Sage script tries for rank n-1 at degree D for TIME_BUDGET_SEC (default 5400s = 1h30m).
#   - If not found, it ALWAYS writes the best-so-far system (possibly deg<D or rank<n-1),
#     plus the classic genlog + secret log.

N=${1:-24}
D=${2:-48}
SEED=${3:--}              # "-" → random (unset)
TIME_BUDGET=${4:-5400}    # 90 minutes internal time budget
MAX_MAPS=${5:-60}         # a bit higher than 40 to use the 90 minutes well
SECRET_TRIES=${6:-4096}   # per-map secret probes
RANK_MIN=${7:--}          # "-" → auto = n-1
FORCE=${8:-false}         # true to overwrite existing .in

mkdir -p logs data/hfe_instances

echo "[$(date)] Node(s): ${SLURM_JOB_NODELIST:-unknown}"
echo "[$(date)] HFE($N,$D) • time_budget=${TIME_BUDGET}s • maps=${MAX_MAPS} • tries/map=${SECRET_TRIES} • rank_min=${RANK_MIN} • force=${FORCE}"

# Let the Sage script choose the default output path by passing "-"
sage scripts/construct_HFE_fast.sage \
  "$N" "$D" "-" "$SEED" "$TIME_BUDGET" "$MAX_MAPS" "$SECRET_TRIES" "$RANK_MIN" "$FORCE"

echo "[$(date)] ✅ Done"
