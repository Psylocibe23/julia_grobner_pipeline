#!/bin/bash
#SBATCH --job-name=solve
#SBATCH --output=logs/%x_%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00
#SBATCH --mem=2G

sage scripts/extract_solutions_from_lex.sage results/your_system_F4_<timestamp>_LEX.txt
