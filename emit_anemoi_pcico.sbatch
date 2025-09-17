#!/usr/bin/env bash
#SBATCH --job-name=emit_pcico
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=logs/slurm/%x_%j.out
#SBATCH --error=logs/slurm/%x_%j.err
# (optional) #SBATCH --partition=cpu

set -euo pipefail

# --- parameters (override via: sbatch --export=ALL,P=251,ALPHA=3,ROUNDS=2 ...) ---
: "${P:=251}" # prime p
: "${ALPHA:=3}" # S-box alpha in {3,5,7,9,11} with gcd(alpha, p-1)=1
: "${ROUNDS:=2}" # number of rounds
: "${REPO_ROOT:=$SLURM_SUBMIT_DIR}" 

cd "$REPO_ROOT"


if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
  source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "$HOME/mambaforge/etc/profile.d/conda.sh" ]; then
  source "$HOME/mambaforge/etc/profile.d/conda.sh"
fi
conda activate gb_conda_env

mkdir -p data logs/slurm

echo "[emit] p=$P alpha=$ALPHA rounds=$ROUNDS"
srun -n 1 sage scripts/emit_anemoi_pcico.sage --p "$P" --alpha "$ALPHA" --rounds "$ROUNDS"

STEM="ANEMOI_p${P}_r${ROUNDS}_a${ALPHA}_PCICO"
echo "[emit] wrote: data/${STEM}.in"
