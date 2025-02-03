#Temporal autocorrelation in year random effects
#Ana Miller-ter Kuile
# June 10, 2024

#this script looks for evidence of temporal auto
#correlation in year random effects for the MSAM and MSOM models

# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 

package.list <- c("here", "tidyverse")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

## And loading them

for(i in package.list){library(i, character.only = T)}


# Function ----------------------------------------------------------------

auto_function <- function(speciesID){
  
  df <- years %>%
    filter(species == speciesID)
  
  auto <- acf(df$`50%`, type = "correlation")
  
  return(auto)
}

#so that more than one plot shows at once
par(mfrow = c(2,2))

# Load data ---------------------------------------------------------------

#grab posteriors for he year random effects for a couple
#of common species and look at whether there is temporal
#correlation in the year effects


# Fish --------------------------------------------------------------------


fish_mod <- readRDS(here('01_sbc_fish',
                    'monsoon',
                    'fish_MSAM',
                    'outputs',
                    'fish_detection_summary.RDS'))

years <- as.data.frame(fish_mod$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(str_detect(parm, "eps.year")) %>%
  filter(parm != 'sig.eps.year') %>%
  filter(str_detect(parm, 'star')) %>%
  separate(parm, 
           into = c('year', 'species'),
           sep = ",") %>%
  mutate(year = str_sub(year, 15, nchar(year))) %>%
  mutate(species = str_sub(species, 1, (nchar(species)-1))) %>%
  mutate(species = as.numeric(species),
         year = as.numeric(year))

species <- unique(years$species)
species_list <- lapply(species, auto_function)

#1111111111111111111111
#22/63 have some temporal correlation 35%

# Birds -------------------------------------------------------------------

bird_mod <- readRDS(here('02_konza_birds',
                         'monsoon',
                         'MSAM',
                         'outputs',
                         'bird_detection_summary.RDS'))

years <- as.data.frame(bird_mod$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(str_detect(parm, "eps.year")) %>%
  filter(parm != 'sig.eps.year') %>%
  filter(str_detect(parm, 'star')) %>%
  separate(parm, 
           into = c('year', 'species'),
           sep = ",") %>%
  mutate(year = str_sub(year, 15, nchar(year))) %>%
  mutate(species = str_sub(species, 1, (nchar(species)-1)))%>%
  mutate(species = as.numeric(species),
         year = as.numeric(year))

species <- unique(years$species)
species_list <- lapply(species, auto_function)

#1111111111111111111111111111111111111111111
#43/78 - 55%


# Hoppers -----------------------------------------------------------------

hop_mod <- readRDS(here('03_sev_grasshoppers',
                        'monsoon',
                        'MSAM',
                        'outputs',
                        'grasshopper_detection_summary.RDS'))

years <- as.data.frame(hop_mod$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(str_detect(parm, "eps.year")) %>%
  filter(parm != 'sig.eps.year') %>%
  filter(str_detect(parm, 'star')) %>%
  separate(parm, 
           into = c('year', 'species'),
           sep = ",") %>%
  mutate(year = str_sub(year, 15, nchar(year))) %>%
  mutate(species = str_sub(species, 1, (nchar(species)-1)))%>%
  mutate(species = as.numeric(species),
         year = as.numeric(year))

species <- unique(years$species)
species_list <- lapply(species, auto_function)

#11111111111111
#14/46 - 30%


# Plants ------------------------------------------------------------------

plant_mod <- readRDS(here('04_nps_plants',
                        'monsoon',
                        'nps_MSAM',
                        'outputs_yrsite',
                        'nps_detection_summary.RDS'))

years <- as.data.frame(plant_mod$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(str_detect(parm, "eps.year")) %>%
  filter(parm != 'sig.eps.year') %>%
  filter(str_detect(parm, 'star')) %>%
  separate(parm, 
           into = c('year', 'species'),
           sep = ",") %>%
  mutate(year = str_sub(year, 15, nchar(year))) %>%
  mutate(species = str_sub(species, 1, (nchar(species)-1)))%>%
  mutate(species = as.numeric(species),
         year = as.numeric(year))

species <- unique(years$species)
species_list <- lapply(species, auto_function)

#1
#1/84 - 1%


