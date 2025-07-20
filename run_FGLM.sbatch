#!/bin/bash
#SBATCH --job-name=fglm
#SBATCH --output=logs/%x_%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem=4G

sage scripts/convert_to_lex_fglm.sage results/your_system_F4_<timestamp>.txt
