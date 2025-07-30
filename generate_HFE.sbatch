#!/bin/bash
#SBATCH --job-name=genhfe           
#SBATCH --output=logs/genhfe_%j.out 
#SBATCH --error=logs/genhfe_%j.err  
#SBATCH --ntasks=1                  
#SBATCH --time=01:00:00             
#SBATCH --mem=8G                    

module load sage

sage scripts/construct_HFE_system.sage
