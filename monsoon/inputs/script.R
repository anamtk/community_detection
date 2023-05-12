#Monsoon script - fish MSOM
# Ana Miller-ter Kuile
# May 12, 2023

#this script runs the fish MSOM

# Load packages ---------------------------------------------------------------
Sys.time()


# Load packages
package.list <- c("jagsUI", "coda") 


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Load Data ---------------------------------------------------------------

#load the formatted data for the JAGS model
data <- readRDS("/scratch/atm234/sbc_fish/inputs/fish_data_dynmultisite.RDS")

# Compile data ------------------------------------------------------------
data_list <- list(y = data$y,
                  z = data$z,
                  vis = data$vis,
                  size = data$size,
                  n.species = data$n.species,
                  n.years = data$n.years,
                  n.start = data$n.start,
                  n.end = data$n.end,
                  n.transects = data$n.transects,
                  n.rep = data$n.rep)

# Parameters to save ------------------------------------------------------

params <- c(
            #COMMUNITY parameters
            "a1.Vis",
            "a2.Size",
            "psi.mean",
            "sd.lpsi",
            "phi.mean",
            "sd.phi",
            'gamma.mean',
            'sd.lgamma',
            'a0.mean',
            'sd.a0')

# INits -------------------------------------------------------------------

#inits <- readRDS("/scratch/atm234/whwo_ipm/parameter_models/adult_occupancy/inputs/adult_occupancy_inits.RDS")

# JAGS model --------------------------------------------------------------

mod <- jagsUI::jags(data = data_list,
                        #inits = inits,
                        inits = NULL,
                        model.file = '/scratch/atm234/sbc_fish/inputs/dynamic_MSOM_multisite.R',
                        parameters.to.save = params,
                        parallel = TRUE,
                        n.chains = 3,
                        n.iter = 4000,
                        DIC = TRUE)

#save as an R data object
saveRDS(mod, 
        file ="/scratch/atm234/sbc_fish/outputs/fish_MSOM_model.RDS")

Sys.time()


