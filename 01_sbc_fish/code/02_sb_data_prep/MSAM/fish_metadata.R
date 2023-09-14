# metadata
# An Bui
# April 2023


# Load packages -----------------------------------------------------------


package.list <- c("here", "tidyverse", "janitor")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# read csv ----------------------------------------------------------------

fish_species <- read_csv(here::here('sbc_fish',
                                    "data_raw", 
                                    "SBCLTER_Species_List_Master_20210113.csv")) %>% 
  clean_names() %>% 
  filter(group == "FISH")

#Export that matrix to a central location for all matrices
saveRDS(fish_species, here('sbc_fish',
                           "data_outputs",
                     "metadata",
                     "species_metadata.RDS"))

