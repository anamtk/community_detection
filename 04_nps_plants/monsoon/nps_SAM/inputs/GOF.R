#Get diss out of converged model
#December 4, 2023

#this script pulls out and summarises diss for visualizations
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

params <- c("diss.rep",
            "resid")

mod2 <- update(mod, 
               parameters.to.save = params,
               n.iter = 1335)


sum2 <- summary(mod2$samples)

saveRDS(sum2, 
        file = "/scratch/sml665/nps_plants/SAM/outputs/nps_SAM_GOF_summary.RDS")
