#Get Bray-Curtis out of converged model
#Ana Miller-ter Kuile
#September 19, 2023

#this script pulls out and sumamrises Bray Curtis for visualizations
#purposes for An's project


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

mod <- readRDS( file ="/scratch/atm234/sbc_fish/outputs/fish_MSAM_model.RDS")

# Update model to track p -------------------------------------------------

parms <- c("p0")

mod2 <- update(mod,
               parameters.to.save = parms,
               n.iter = 4000)


# Summary of spp level detectoin ------------------------------------------


sum2 <- summary(mod2$samples)

saveRDS(sum2, 
        file = "/scratch/atm234/sbc_fish/outputs/fish_p0_summary.RDS")


# Get summary of OG model -------------------------------------------------

sum <- summary(mod$samples)

saveRDS(sum, 
        file = "/scratch/atm234/sbc_fish/outputs/fish_detection_summary.RDS")
