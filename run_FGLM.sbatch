#!/bin/bash
#SBATCH --job-name=fglm
#SBATCH --output=logs/%x_%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem=32G

# Usage: sbatch fglm.sbatch results/HFE_n4_d4_system_F4_20250731_000946.txt

INPUT_FILE="$1"
if [[ -z "$INPUT_FILE" ]]; then
    echo "Usage: sbatch fglm.sbatch <F4_result_file>"
    exit 1
fi

module load sagemath/10.6-singularity
sage -python scripts/convert_to_lex_fglm.sage.py "$INPUT_FILE"
