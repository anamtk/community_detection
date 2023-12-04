#Get diss out of converged model
#December 4, 2023

#this script pulls out and sumamrises Bray Curtis for visualizations
#purposes for An's project


# Load packages ---------------------------------------------------------------
Sys.time()


# Load packages,
package.list <- c('jagsUI',"coda") 


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load model --------------------------------------------------------------

mod <- readRDS( file = "/scratch/sml665/nps_plants/SAM/outputs/nps_SAM_model.RDS")


# Get summary of model -------------------------------------------------

sum <- summary(mod$samples)

saveRDS(sum, 
        file = "/scratch/sml665/nps_plants/SAM/outputs/nps_SAM_summary.RDS")


