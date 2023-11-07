#Convergence check for nps MSAM
#Shelby Lamm
#October 30, 2023

#this script assesses convergence in the nps msam

# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 
package.list <- c("here", "tidyverse") 


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

source(here('00_functions',
            'plot_functions.R'))

# Load Data ---------------------------------------------------------------

rhat <- readRDS(here('04_nps_plants',
                     'monsoon',
                     "nps_MSAM",
                     'outputs',
                     'nps_MSAM_model_Rhat4.RDS'))

rhat <- readRDS(here('04_nps_plants',
                     'monsoon',
                     "nps_MSAM",
                     'outputs_no_lifegroup',
                     'nps_MSAM_model_Rhat3.RDS'))

# Graph RHat per parameter ------------------------------------------------

rhat_graph_fun(rhat)
