#!/bin/bash
#SBATCH --job-name=gb_F4
#SBATCH --output=logs/%x_%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --mem=8G

julia scripts/solve_F4_from_file.jl data/your_system.in 8
