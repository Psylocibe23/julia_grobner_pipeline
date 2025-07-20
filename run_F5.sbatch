#!/bin/bash
#SBATCH --job-name=gb_F5
#SBATCH --output=logs/%x_%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem=8G

julia scripts/solve_F5_from_file.jl data/your_system.in
