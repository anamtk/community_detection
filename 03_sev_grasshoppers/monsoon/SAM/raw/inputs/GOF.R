#Get Bray-Curtis out of converged model
#Ana Miller-ter Kuile
#September 19, 2023

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

mod <- readRDS( file = "/scratch/atm234/sev_hoppers/SAM/raw/outputs/sev_SAM_model_raw.RDS")


# Get summary of model -------------------------------------------------

params <- c("bray.rep",
            "resid")

mod2 <- update(mod, 
               parameters.to.save = params,
               n.iter = 1335)


sum2 <- summary(mod2$samples)

saveRDS(sum2, 
        file = "/scratch/atm234/sev_hoppers/SAM/raw/outputs/hopper_SAM_raw_GOF_summary.RDS")
