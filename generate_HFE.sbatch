#!/bin/bash
#SBATCH --job-name=genhfe           
#SBATCH --output=logs/genhfe_%j.out 
#SBATCH --error=logs/genhfe_%j.err  
#SBATCH --partition=turing-long
#SBATCH --ntasks=1                  
#SBATCH --time=02:00:00             
#SBATCH --cpus-per-task=32           
#SBATCH --mem=128G                                 


module load sagemath/10.6-singularity

# Run the generator script with Sage
echo "Starting HFE(n,d) generation on $SLURM_JOB_NODELIST"
sage -python scripts/construct_HFE_system.sage.py
echo "✅ Done at $(date)"
