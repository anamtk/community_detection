# Running the nest initiation
# Ana Miller-ter Kuile
# November 4, 2021

# Load packages ---------------------------------------------------------------
Sys.time()
# Load packages, here and tidyverse for coding ease, 
package.list <- c("jagsUI") #mcmc output


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Load Data ---------------------------------------------------------------


data_list <- readRDS("/scratch/atm234/sbc_benthic/inputs/input_data.RDS")
data <- readRDS(here('monsoon', 'inputs', 'input_data.RDS'))
data <- list(y = data_list$y,
             n.species = data_list$n.species,
             n.years = data_list$n.years,
             n.reps = data_list$n.reps,
             vis = data_list$vis,
             z = data_list$z)
# Parameters to save ------------------------------------------------------

params <- c(#species-level parameters
            'psi',
            "eta",
            "a1.Vis",
            #community-level parameters
            'tau.eta',
            'sd.vis',
            'sd.lpsi',
            'sd.lp',
            "rho"
            )


# JAGS model --------------------------------------------------------------

jags <- jagsUI::jags(data = data,
                     inits = NULL,
                     model.file = '/scratch/atm234/sbc_benthic/inputs/model.R',
                     parameters.to.save = params,
                     parallel = TRUE,
                     n.chains = 3,
                     n.iter = 70000,
                     n.burnin = 1000,
                     DIC = TRUE)

saveRDS(jags,
        file = "/scratch/atm234/sbc_benthic/outputs/benthic_MSOM_1_20.RDS")

Sys.time()




