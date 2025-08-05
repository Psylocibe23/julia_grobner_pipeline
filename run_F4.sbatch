#!/bin/bash
#SBATCH --job-name=gb_F4
#SBATCH --output=logs/%x_%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --mem=32G

# Usage: sbatch gb_F4.sbatch data/HFE_n4_d4_system.in 8

INPUT_FILE="$1"   # e.g., data/HFE_n4_d4_system.in
THREADS="${2:-8}" # Defaults to 8 if not given

if [[ -z "$INPUT_FILE" ]]; then
    echo "Usage: sbatch gb_F4.sbatch <input_file> [threads]"
    exit 1
fi

module load julia
julia scripts/solve_F4_from_file.jl "$INPUT_FILE" "$THREADS"
