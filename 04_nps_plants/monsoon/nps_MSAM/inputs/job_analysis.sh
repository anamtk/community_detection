#!/bin/bash
#SBATCH --job-name=nps_MSAM_subset_analysis
#SBATCH --output=/scratch/sml665/nps_plants/outputs/modelmsam_analysis.out
#SBATCH --time=3:00:00                                                               
#SBATCH --chdir=/scratch/sml665/nps_plants/inputs
#SBATCH --mem=85000                                                                    
#SBATCH --mail-type=all
#SBATCH --mail-user=sml665@nau.edu
#SBATCH --cpus-per-task=1


module load R/4.1.2
module  --ignore-cache load jags/4.3.0

srun Rscript MSAM_script_analysis.R