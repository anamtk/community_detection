
# Load packages ---------------------------------------------------------------
Sys.time()

# Load packages,
#package.list <- c('jagsUI',"coda",'dplyr', 
#                  'stringr','magrittr',
#                  'tidyr','ggplot2',
#                  'tibble', 'reshape2',
#                 'data.table') 

library(jagsUI)
library(coda)
library(dplyr)
library(purrr)
library(stringr)
library(magrittr)
library(tidyr)
library(ggplot2)
library(tibble)
library(reshape2)
library(data.table)
## Installing them if they aren't already on the computer
#new.packages <- package.list[!(package.list %in% 
#                                 installed.packages()[,"Package"])]
#if(length(new.packages)) install.packages(new.packages)

## And loading them
#for(i in package.list){library(i, character.only = T)}



# Load model --------------------------------------------------------------

model <- readRDS(file ="/scratch/atm234/nps_plants/outputs/nps_JAGS_RE_model.RDS")

print('model loaded')

# update for z ------------------------------------------------------------

Sys.time()
print("start z model")

parms2 <- c("z")

z_update <- update(model, 
                   parameters.to.save = parms2,
                   n.iter = 700)

Sys.time()
print("end z model")

sum <- summary(s_update$samples)

saveRDS(sum, file ="/scratch/atm234/nps_plants/outputs/plant_msom_z_sum.RDS")