#Get GOF of MSAM
#Ana Miller-ter Kuile
#September 19, 2023

#this script pulls out and sumamrises y.rep and N for the MSAM


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

mod <- readRDS( file = "/scratch/atm234/sev_hoppers/outputs/sev_MSAM_nocovRE_model.RDS")


# Get summary of model -------------------------------------------------

params <- c("y.rep",
            "N")

mod2 <- update(mod, 
               parameters.to.save = params,
               n.iter = 335)


sum2 <- summary(mod2$samples)

saveRDS(sum2, 
        file = "/scratch/atm234/sev_hoppers/outputs/sev_MSAM_GOF_summary.RDS")
