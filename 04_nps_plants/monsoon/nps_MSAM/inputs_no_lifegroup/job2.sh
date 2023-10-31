#!/bin/bash
#SBATCH --job-name=nps_MSAM_subset_no_lifegroup2
#SBATCH --output=/scratch/sml665/nps_plants/outputs_no_lifegroup/modelmsam_subset2.out
#SBATCH --time=24:00:00                                               
#SBATCH --chdir=/scratch/sml665/nps_plants/inputs_no_lifegroup
#SBATCH --mem=85000                                                                    
#SBATCH --mail-type=all
#SBATCH --mail-user=sml665@nau.edu
#SBATCH --cpus-per-task=3


module load R/4.1.2
module  --ignore-cache load jags/4.3.0

srun Rscript MSAM_script2.R