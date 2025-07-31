#!/bin/bash
#SBATCH --job-name=fglm
#SBATCH --output=logs/%x_%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem=4G

sage scripts/convert_to_lex_fglm.sage results/HFE_n4_d4_system_F4_20250731_000946.txt
