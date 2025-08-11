#!/bin/bash
#SBATCH --job-name=gen_hfe
#SBATCH --output=logs/%x_%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --mem=64G

# Usage:
#   sbatch run_generate_HFE.sbatch <n> <D> [outfile.in] [seed] [prob_quad] [prob_lin] [rank_min] [allow_fallback] [max_maps] [max_secret_tries]
# Example:
#   sbatch run_generate_HFE.sbatch 80 96 - 123 0.7 0.6 79 true 30 4096

set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: sbatch run_generate_HFE.sbatch <n> <D> [outfile.in] [seed] [prob_quad] [prob_lin] [rank_min] [allow_fallback] [max_maps] [max_secret_tries]"
  exit 1
fi

module load sagemath/10.6-singularity
sage scripts/construct_HFE_system.sage "$@"
