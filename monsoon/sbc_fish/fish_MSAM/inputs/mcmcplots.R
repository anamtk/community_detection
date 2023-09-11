#Monsoon script - fish MSAM
# Ana Miller-ter Kuile
# July 27, 2023

#this script runs the fish MSOM

# Load packages ---------------------------------------------------------------
Sys.time()


# Load packages,
package.list <- c("jagsUI", "coda",
                  'dplyr', 'stringr',
                  'magrittr', 'tidyr',
                  'mcmcplots','ggplot2') 


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Load Data ---------------------------------------------------------------

#load the formatted data for the JAGS model
model <- readRDS(file ="/scratch/atm234/sbc_fish/outputs/fish_MSAM_model.RDS")

params <- c(
  #COMMUNITY parameters
  'a1.Vis',
  'a2.Size',
  'lambda.mean',
  'sig.lambda',
  'a0.mean',
  'sig.a0')

# Check convergence -------------------------------------------------------

mcmcplot(model$samples,
         dir = "/scratch/atm234/sbc_fish/outputs/mcmcplots/MSAM",
         parms = params)


