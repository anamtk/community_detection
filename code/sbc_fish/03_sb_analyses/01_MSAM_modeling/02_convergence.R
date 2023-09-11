#Convergence check for fish MSAM
# Ana Miller-ter Kuile
#Sepptember 7, 2023

#this script assesses convergence in the fish msam

# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 
package.list <- c("here", "tidyverse") 


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load Data ---------------------------------------------------------------

rhat <- readRDS(here('monsoon',
                     'sbc_fish',
                     'fish_MSAM',
                     'outputs',
                     'fish_MSAM_model_Rhat.RDS'))


rhat2 <- readRDS(here('monsoon',
                     'sbc_fish',
                     'fish_MSAM',
                     'outputs',
                     'fish_MSAM_model_Rhat_update.RDS'))

# Graph RHat per parameter ------------------------------------------------

rhat_graph_fun(rhat)
rhat_graph_fun(rhat2)
