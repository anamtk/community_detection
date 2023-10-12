#!/bin/bash
#SBATCH --job-name=nps_MSAM
#SBATCH --output=/scratch/sml665/nps_plants/outputs/modelmsam.out
#SBATCH --time=48:00:00                                                               
#SBATCH --chdir=/scratch/sml665/nps_plants/inputs
#SBATCH --mem=85000                                                                    
#SBATCH --mail-type=all
#SBATCH --mail-user=sml665@nau.edu


module load R/4.1.2
module load jags/4.3.0

srun Rscript MSAM_script.R