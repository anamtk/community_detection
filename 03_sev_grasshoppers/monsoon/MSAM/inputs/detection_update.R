#Get baseline detection out of converged model
#Ana Miller-ter Kuile
#September 19, 2023

#this script pulls out and sumamrises detection for visualizations
#and potential further analyses

# Load packages ---------------------------------------------------------------
Sys.time()


# Load packages,
package.list <- c('jagsUI',"coda",'dplyr', 
                   'stringr','magrittr',
                   'tidyr','ggplot2',
                   'tibble') 


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load model --------------------------------------------------------------

mod <- readRDS( file ="/scratch/atm234/sev_hoppers/outputs/sevMSAMnocov_model2.RDS")

# Update model to track p -------------------------------------------------

parms <- c("p")

mod2 <- update(mod,
               parameters.to.save = parms,
               n.iter = 4000)


# Summary of spp level detectoin ------------------------------------------


sum2 <- summary(mod2$samples)

saveRDS(sum2, 
        file = "/scratch/atm234/sev_hoppers/outputs/grasshopper_p_summary.RDS")


