#!/bin/bash
#SBATCH --job-name=genhfe           
#SBATCH --output=logs/genhfe_%j.out 
#SBATCH --error=logs/genhfe_%j.err  
#SBATCH --partition=turing-long
#SBATCH --ntasks=1                  
#SBATCH --time=08:00:00             
#SBATCH --cpus-per-task=32          
#SBATCH --mem=128G                                 

module load sagemath/10.6-singularity

# Example values for n and d
N=6
D=8

echo "Starting HFE($N,$D) generation on $SLURM_JOB_NODELIST"
sage scripts/construct_HFE_parallelized.sage $N $D
echo "✅ Done at $(date)"