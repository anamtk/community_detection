#Convergence check for MSOM
# Ana Miller-ter Kuile
#September 7, 2023

#this script assesses convergence in an MSOM/MSAM downloaded from a monsoon run

# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 
package.list <- c("here", "tidyverse") 


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

source(here("00_functions",
            'plot_functions.R'))

# Load Data ---------------------------------------------------------------

#this is the rhat that you can download with the MSAM_Script.R from monsoon
#eg this file: "/scratch/atm234/konza_birds/outputs/bird_MSAM_model_Rhat.RDS"

rhat <- readRDS(here('04_nps_plants',
                     'example',
                     'monsoon',
                     'outputs',
                     'bird_MSAM_model_Rhat.RDS'))


# Graph RHat per parameter ------------------------------------------------

rhat_graph_fun(rhat)

