#!/bin/bash
#SBATCH --job-name=nps_SAM_subset
#SBATCH --output=/scratch/sml665/nps_plants/SAM/outputs/modelsam_subset.out
#SBATCH --time=24:00:00                                                               
#SBATCH --chdir=/scratch/sml665/nps_plants/SAM/inputs
#SBATCH --mem=85000                                                                    
#SBATCH --mail-type=all
#SBATCH --mail-user=sml665@nau.edu
#SBATCH --cpus-per-task=3


module load R/4.1.2
module  --ignore-cache load jags/4.3.0

srun Rscript SAM_script.R