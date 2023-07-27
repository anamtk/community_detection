#Monsoon script - fish MSAM
# Ana Miller-ter Kuile
# July 27, 2023

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
data <- readRDS("/scratch/atm234/sbc_fish/inputs/fish_msam_dynmultisite.RDS")

# Compile data ------------------------------------------------------------
data_list <- list(y = data$y,
             vis = data$vis,
             size = data$size,
             n.species = data$n.species,
             n.years = data$n.years,
             n.start = data$n.start,
             n.end = data$n.end,
             n.transects = data$n.transects,
             n.rep = data$n.rep,
             #for initials
             ymax = data$ymax)

# Parameters to save ------------------------------------------------------

params <- c(
  #COMMUNITY parameters
  'a1.Vis',
  'a2.Size',
  'lambda.mean',
  'sig.lambda',
  'a0.mean',
  'sig.a0')

# INits -------------------------------------------------------------------

#we found ymax to set initials, since otherwise the model will hate us
inits <- function() list(N = data$ymax)

# JAGS model --------------------------------------------------------------

mod <- jagsUI::jags(data = data_list,
                        inits = inits,
                        #inits = NULL,
                        model.file = '/scratch/atm234/sbc_fish/inputs/MSAM_multisite.R',
                        parameters.to.save = params,
                        parallel = TRUE,
                        n.chains = 3,
                        n.iter = 20000,
                        n.burnin = 100,
                        DIC = TRUE)

#save as an R data object
saveRDS(mod, 
        file ="/scratch/atm234/sbc_fish/outputs/fish_MSAM_model.RDS")

Sys.time()



# Check convergence -------------------------------------------------------

mcmcplot(mod$samples,
         dir = "/scratch/atm234/sbc_fish/outputs/mcmcplots/MSAM")

# Get RHat per parameter ------------------------------------------------

Rhat <- mod$Rhat

saveRDS(Rhat, "/scratch/atm234/sbc_fish/outputs/fish_MSAM_model_Rhat.RDS")



