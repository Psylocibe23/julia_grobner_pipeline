#!/bin/bash
#SBATCH --job-name=genhfe           
#SBATCH --output=logs/genhfe_%j.out 
#SBATCH --error=logs/genhfe_%j.err  
#SBATCH --ntasks=1                  
#SBATCH --time=02:00:00             
#SBATCH --mem=32G                    

module load sagemath/10.6-singularity


sage scripts/construct_HFE_system.sage
