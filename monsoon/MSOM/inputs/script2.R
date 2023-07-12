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


# Load model --------------------------------------------------------------


#load the formatted data for the JAGS model
mod <- readRDS("/scratch/atm234/sbc_fish/outputs/fish_MSOM_model.RDS")

# Parameters to save ------------------------------------------------------

params <- c(
            #COMMUNITY parameters
            "A",
            "B",
            "C")

# INits -------------------------------------------------------------------

mod2 <- update(mod,
               parameters.to.save = params,
               n.iter = 500)


#save as an R data object
saveRDS(mod2, 
        file ="/scratch/atm234/sbc_fish/outputs/fish_MSOM_model_2.RDS")

Sys.time()


