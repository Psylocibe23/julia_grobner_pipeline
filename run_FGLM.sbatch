#!/bin/bash
#SBATCH --job-name=fglm
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --partition=turing-long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --chdir=/home/fspreafichi/julia_grobner_pipeline

set -euo pipefail
set -x

if [[ $# -lt 1 ]]; then
  echo "Usage: sbatch run_FGLM.sbatch <F4_result_file.txt>"
  exit 1
fi
INPUT_FILE="$1"
if [[ ! -f "$INPUT_FILE" ]]; then
  echo "ERROR: file not found: $INPUT_FILE"
  exit 2
fi

module load sagemath/10.6-singularity

# Avoid BLAS/OMP oversubscription (Singular is single-threaded here)
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# For Boolean+field-eq systems we can assume 0-dim safely.
sage scripts/convert_to_lex_fglm.sage "$INPUT_FILE" \
  --assume-zerodim \
  --fglm-timeout=3600 \
  --reduce=auto \
  --reduce-timeout=300
