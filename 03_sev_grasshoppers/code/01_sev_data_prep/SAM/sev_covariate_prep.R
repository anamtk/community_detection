#Sevilleta grasshoppers environmental variable prep
#Ana Miller-ter Kuile
# October 12, 2023

#this script preps the environmental variables with lags
#for the Sevilleta Grasshopper dataset


#Sites for Sevilleta grasshoppers, based on temporal coverage:

#NPP is from 1999 for two sites: 
#Core-Blak gramma and Core-Creosote
#these sites are also called "Five Points" in the datasets
#I think for meteorological station this includes station:
#49 (Five points) 

#NPP are seasonal - spring and fall of each year
#looks like grouping them by "web" and then taking an average of 
#the plots on that web for all species biomass in each season
#would be the way to combine with grasshopper data

#Climate are hourly - only one site near both the creosote
#and black gramma sites
#I think for meteorological station this includes station:
#49 (Five points) 

# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 'data.table')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Load data ---------------------------------------------------------------

c1 <- read.csv(here('03_sev_grasshoppers',
                    'data_raw',
                    'environmental',
                    'Sevilleta_LTER_Hourly_Meteorological_Data_1995_1999.csv'))

c2 <- read.csv(here('03_sev_grasshoppers',
                    'data_raw',
                    'environmental',
                    'Sevilleta_LTER_Hourly_Meteorological_Data_2000_2004.csv'))

c3 <- read.csv(here('03_sev_grasshoppers',
                    'data_raw',
                    'environmental',
                    'Sevilleta_LTER_Hourly_Meteorological_Data_2005_2009.csv'))

c4 <- read.csv(here('03_sev_grasshoppers',
                    'data_raw',
                    'environmental',
                    'Sevilleta_LTER_Hourly_Meteorological_Data_2010_2014.csv'))

c5 <- read.csv(here('03_sev_grasshoppers',
                    'data_raw',
                    'environmental',
                    'Sevilleta_LTER_Hourly_Meteorological_Data_2015_2019.csv'))

c6 <- read.csv(here('03_sev_grasshoppers',
                    'data_raw',
                    'environmental',
                    'Sevilleta_LTER_Hourly_Meteorological_Data_2020_2022.csv'))

npp <- read.csv(here('03_sev_grasshoppers',
                     'data_raw',
                     'environmental',
                     'sev331_quadrat_plant_species_biomass.csv'))

# Combine all climate datasets --------------------------------------------

#combine and select only the weather station we want
climate <- rbind(c1, c2, c3, c4, c5, c6) %>%
  filter(StationID == 49)

temp <- climate %>%
  group_by(Month, Year) %>%
  summarise(Temp = mean(Temp_C, na.rm = T)) %>%
  ungroup()

ppt <- climate %>%
  group_by(Month, Year) %>%
  summarise(PPT = mean(Precipitation, na.rm = T)) %>%
  ungroup()


# NPP ---------------------------------------------------------------------

#NPP are seasonal - spring and fall of each year
#looks like grouping them by "web" and then taking an average of 
#the plots on that web for all species biomass in each season
#would be the way to combine with grasshopper data

#could also look at sites with < 20 years, which would be 
#core_blue and core_PJ, which would be 
#met stations 50 and 42, if we add them in would need
#to also include those weather data

npp %>%
  group_by(site) %>%
  slice_min(order_by = year) %>%
  distinct(site, year, MetStation) %>%
  arrange(year)

npp2 <- npp %>%
  group_by(site, year, season,web, plot, quad ) %>%
  summarise(biomass = sum(biomass.BM, na.rm = T)) %>%
  ungroup() %>%
  filter(site %in% c("core_black", "core_blue"))

npp3 <- npp2 %>%
  group_by(site, web, year, season) %>%
  summarise(NPP = mean(biomass, na.rm = T)) %>%
  ungroup() %>%
  mutate(seasonnum = case_when(season == "spring" ~ 1,
                               season == "fall" ~ 2,
                               TRUE ~ NA_real_))

# Get lags ----------------------------------------------------------------

temp_lags <- temp %>%
  arrange(Year, Month) %>%
  #this creates a column for every lag this month to 11 months ago
  do(data.frame(., setNames(shift(.$Temp, 1:11), c("Temp_l1", 'Temp_l2', "Temp_l3",
                                                 "Temp_l4", "Temp_l5",
                                                 'Temp_l6', 'Temp_l7',
                                                 'Temp_l8', 'Temp_l9',
                                                 'Temp_l10', "Temp_l11")))) %>%
  ungroup()

ppt_lags <- ppt %>%
  arrange(Year, Month) %>%
  #this creates a column for every lag this month to 11 months ago
  do(data.frame(., setNames(shift(.$PPT, 1:11), c("PPT_l1", 'PPT_l2', "PPT_l3",
                                                   "PPT_l4", "PPT_l5",
                                                   'PPT_l6', 'PPT_l7',
                                                   'PPT_l8', 'PPT_l9',
                                                   'PPT_l10', "PPT_l11")))) %>%
  ungroup()

#get npp lags:
npp_lags <- npp3 %>%
  group_by(site) %>%
  arrange(year, seasonnum)  %>%
  #this creates a column for every lag this year to 10 seasons ago
  do(data.frame(., setNames(shift(.$NPP, 1:10), c("NPP_l1", 'NPP_l2', "NPP_l3",
                                                 "NPP_l4", "NPP_l5", "NPP_l6",
                                                 "NPP_l7", "NPP_l8", "NPP_l9",
                                                 "NPP_l10")))) %>%
  ungroup() %>%
  #consider "this season" e.g. fall - and then previous seasons from
  #that for each web
  filter(seasonnum == 2) %>%
  dplyr::select(-seasonnum)
