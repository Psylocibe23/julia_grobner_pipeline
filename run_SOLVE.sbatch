#!/bin/bash
#SBATCH --job-name=solve
#SBATCH --output=logs/%x_%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem=32G

# Usage: sbatch solve.sbatch results/HFE_n4_d4_system_F4_20250731_000946_LEX.txt

INPUT_FILE="$1"
if [[ -z "$INPUT_FILE" ]]; then
    echo "Usage: sbatch solve.sbatch <LEX_basis_file>"
    exit 1
fi

module load sagemath/10.6-singularity
sage scripts/extract_solutions_from_lex.sage "$INPUT_FILE"
