#!/bin/bash
#SBATCH --job-name=fglm_adj
#SBATCH --partition=turing-long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1         
#SBATCH --mem=64G                
#SBATCH --time=08:00:00           
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Usage:
#   sbatch run_FGLM_adjusted.sbatch results/HFE_n40_D80_F4_20250828_020103.txt \
#           --assume-zerodim --fglm-timeout=1800 --reduce=auto --reduce-timeout=300

set -euo pipefail

mkdir -p logs results

if [[ $# -lt 1 ]]; then
  echo "Usage: sbatch run_FGLM_adjusted.sbatch <F4_or_F5_result.txt> [extra flags...]"
  exit 1
fi

INPUT_FILE="$1"; shift || true
EXTRA_ARGS=("$@")

if [[ ! -f "$INPUT_FILE" ]]; then
  echo "ERROR: input file not found: $INPUT_FILE"
  exit 2
fi

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "[$(date)] Node(s): ${SLURM_JOB_NODELIST:-unknown}"
echo "[$(date)] CWD:      $(pwd)"
echo "[$(date)] Input:    $INPUT_FILE"
echo "[$(date)] Flags:    ${EXTRA_ARGS[*]:-<none>}"

# Make 'module' available in non-login shells, then load Sage 
if ! type module >/dev/null 2>&1; then
  [[ -f /etc/profile ]] && source /etc/profile
  [[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh
fi
module load sagemath/10.6-singularity

# Run the adjusted FGLM driver
echo "[$(date)] Launching Sage…"
set -x
sage scripts/fglm_adjusted.sage "$INPUT_FILE" "${EXTRA_ARGS[@]}"
set +x

stem="$(basename "$INPUT_FILE" .txt)"
lex_out="results/${stem}_LEX.txt"
log_out="logs/${stem}_FGLM_adjusted.log"

echo "[$(date)] Expected outputs:"
echo "  LEX basis : ${lex_out}   $( [[ -f "$lex_out" ]] && echo '[OK]' || echo '[missing]' )"
echo "  FGLM log  : ${log_out}   $( [[ -f "$log_out" ]] && echo '[OK]' || echo '[missing]' )"

echo "[$(date)] Listing results/:"
ls -lh results || true

echo "[$(date)] ✅ Done"
