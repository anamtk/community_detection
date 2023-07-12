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
data <- readRDS("/scratch/atm234/sbc_fish/SAM/inputs/turnover_SAM_input_data.RDS")

# Compile data ------------------------------------------------------------
data_list <- list(n.data = data$n.data,
                  n.transects = data$n.transects,
                  n.sites = data$n.sites,
                  n.years = data$n.years,
                  Transect.ID = data$Transect.ID,
                  Site.ID = data$Site.ID,
                  Year.ID = data$Year.ID,
                  gain = data$gain,
                  n.kelplag = data$n.kelplag,
                  Kelp = data$Kelp,
                  n.templag = data$n.templag,
                  Temp = data$Temp,
                  Chla = data$Chla)

# Parameters to save ------------------------------------------------------

params <- c("b0",
            'b',
            'wA',
            'wB',
            'wC')

# INits -------------------------------------------------------------------

#inits <- readRDS("/scratch/atm234/whwo_ipm/parameter_models/adult_occupancy/inputs/adult_occupancy_inits.RDS")

# JAGS model --------------------------------------------------------------

mod <- jagsUI::jags(data = data_list,
                        #inits = inits,
                        inits = NULL,
                        model.file = '/scratch/atm234/sbc_fish/SAM/inputs/gain_SAM.R',
                        parameters.to.save = params,
                        parallel = TRUE,
                        n.chains = 3,
                        n.iter =  51027,
                        n.burnin = 1000,
                        DIC = TRUE)

#save as an R data object
saveRDS(mod, 
        file ="/scratch/atm234/sbc_fish/SAM/outputs/gain_model.RDS")

Sys.time()


