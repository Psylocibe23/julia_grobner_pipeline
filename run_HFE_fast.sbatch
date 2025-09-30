#!/bin/bash
#SBATCH --job-name=hfe_fast
#SBATCH --output=logs/hfe_fast_%j.out
#SBATCH --error=logs/hfe_fast_%j.err
#SBATCH --qos=tesi_spreafichi
#SBATCH --partition=turing-wide
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8       # keep modest to avoid QOSGrpCpuLimit
#SBATCH --mem=32G               # sufficient for generation
#SBATCH --time=02:00:00

set -euo pipefail

module load sagemath/10.6-singularity

# Usage:
#   sbatch run_HFE_fast.sbatch N D [SEED] [TIME_BUDGET_SEC] [MAX_MAPS] [SECRET_TRIES] [RANK_MIN|-] [FORCE]
#
# Typical:
#   sbatch run_HFE_fast.sbatch 24 48 0 5400 40 2048 - true

N=${1:-24}
D=${2:-48}
SEED=${3:--}
TIME_BUDGET=${4:-5400}        # 90 minutes inside the job
MAX_MAPS=${5:-40}
SECRET_TRIES=${6:-2048}
RANK_MIN=${7:--}              # "-" → auto n-1
FORCE=${8:-false}

mkdir -p logs data/hfe_instances

echo "[$(date)] Node(s): ${SLURM_JOB_NODELIST:-unknown}"
echo "[$(date)] HFE($N,$D) • time_budget=${TIME_BUDGET}s • maps=${MAX_MAPS} • tries/map=${SECRET_TRIES} • rank_min=${RANK_MIN} • force=${FORCE}"

# Let the script choose the default output path (pass "-")
sage scripts/construct_HFE_fast.sage \
     "$N" "$D" "-" "$SEED" "$TIME_BUDGET" "$MAX_MAPS" "$SECRET_TRIES" "$RANK_MIN" "$FORCE"

echo "[$(date)] ✅ Done"
