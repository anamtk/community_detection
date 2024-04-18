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

#to get the sites and transects we need for the other 
#two datasets
stability <- readRDS(here("03_sev_grasshoppers",
                          "monsoon",
                          "MSAM",
                          "outputs",
                          "sev_bray_meanSD.RDS"))

IDs <- read.csv(here('03_sev_grasshoppers',
                     'data_outputs',
                     'metadata',
                     'sev_site_year_IDs.csv'))

raw_dat <- readRDS(here('05_visualizations',
                        'viz_data',
                        'sev_observed_bray.RDS'))

# Get response data in order ----------------------------------------------

stability2 <- as.data.frame(stability) %>%
  rownames_to_column(var = "var") %>%
  filter(var != "deviance") %>%
  separate(var,
           into = c('siteID', 'yrID'),
           sep = ",") %>%
  mutate(siteID = str_sub(siteID, 6, nchar(siteID))) %>%
  mutate(yrID = str_sub(yrID, 1, (nchar(yrID)-1))) %>%
  rename("bray" = "Mean") %>%
  dplyr::select(yrID, siteID, bray, SD) %>%
  mutate(yrID = as.numeric(yrID),
         siteID = as.numeric(siteID)) %>%
  left_join(IDs, by = c("siteID", "yrID"))

# Combine all climate datasets --------------------------------------------

#combine and select only the weather station we want
climate <- rbind(c1, c2, c3, c4, c5, c6) %>%
  filter(StationID == 49)


# Temperature -------------------------------------------------------------


#monthly temps
temp <- climate %>%
  group_by(Month, Year) %>%
  summarise(Temp = mean(Temp_C, na.rm = T)) %>%
  ungroup()

ggplot(temp, aes(x = Month, y= Temp)) +
  geom_jitter(width = 0.2, height = 0)

#which are above and below mean?
climate %>%
  mutate(Temp_C = scale(Temp_C)) %>%
  group_by(Month) %>%
  summarise(Temp = mean(Temp_C, na.rm = T)) %>%
  ungroup()
#1, 2, 3, 11, 12 - cold
#4, 5, 6, 7, 8, 9, 10 - warm

#seasonal temperatures - 1 is "cold", 2 is "warm"
temp2 <- climate %>%
  mutate(season = case_when(Month %in% c(1, 2, 3, 11, 12) ~ 1,
                            Month %in% c(4, 5, 6, 7, 8, 9, 10) ~ 2))  %>%
  group_by(season, Year) %>%
  summarise(Temp = mean(Temp_C, na.rm = T)) %>%
  ungroup()


# Precip ------------------------------------------------------------------

ppt <- climate %>%
  group_by(Month, Year) %>%
  summarise(PPT = mean(Precipitation, na.rm = T)) %>%
  ungroup()

ggplot(ppt, aes(x = Month, y= PPT)) +
  geom_jitter(width = 0.2, height = 0)

climate %>%
  mutate(PPT = scale(Precipitation)) %>%
  group_by(Month) %>%
  summarise(PPT = mean(PPT, na.rm = T)) %>%
  ungroup()

#1, 2, 3, 4, 5, 6, 11, 12 - dry
#7, 8, 9, 10 - wet

#seasonal ppt 1 is dry, 2 is wet
ppt2 <- climate %>%
  mutate(season = case_when(Month %in% c(1, 2, 3, 4, 5, 6, 11, 12) ~ 1,
                            Month %in% c(7, 8, 9, 10) ~ 2))  %>%
  group_by(season, Year) %>%
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
  filter(site %in% c("core_black", "core_creosote"))

npp3 <- npp2 %>%
  group_by(site, web, year, season) %>%
  summarise(NPP = mean(biomass, na.rm = T)) %>%
  ungroup() %>%
  mutate(seasonnum = case_when(season == "spring" ~ 1,
                               season == "fall" ~ 2,
                               TRUE ~ NA_real_))

# Get lags ----------------------------------------------------------------

#remove all months except the last survey month. which either SEptember or OCtober
#depending on year and site. Let's just go ahead and set October
temp_lags <- temp %>%
  arrange(Year, Month) %>%
  #this creates a column for every lag this month to 11 months ago
  do(data.frame(., setNames(shift(.$Temp, 1:11), c("Temp_l1", 'Temp_l2', "Temp_l3",
                                                 "Temp_l4", "Temp_l5",
                                                 'Temp_l6', 'Temp_l7',
                                                 'Temp_l8', 'Temp_l9',
                                                 'Temp_l10', "Temp_l11")))) %>%
  ungroup() %>%
  #start in october, which is when sampling ended
  filter(Month == 10) %>%
  dplyr::select(-Month)

#seasonal lags
temp_lags2 <- temp2 %>%
  arrange(Year, season) %>%
  #this creates a column for every lag this season to 5 ago
  do(data.frame(., setNames(shift(.$Temp, 1:5), c("Temp_l1", 'Temp_l2', "Temp_l3",
                                                   "Temp_l4", "Temp_l5")))) %>%
  ungroup() %>%
  #start in cold, when surveys started
  filter(season == 1) %>%
  dplyr::select(-season)
  

ppt_lags <- ppt %>%
  arrange(Year, Month) %>%
  #this creates a column for every lag this month to 11 months ago
  do(data.frame(., setNames(shift(.$PPT, 1:11), c("PPT_l1", 'PPT_l2', "PPT_l3",
                                                   "PPT_l4", "PPT_l5",
                                                   'PPT_l6', 'PPT_l7',
                                                   'PPT_l8', 'PPT_l9',
                                                   'PPT_l10', "PPT_l11")))) %>%
  ungroup() %>%
  #start in october, which is when sampling ended
  filter(Month == 10) %>%
  dplyr::select(-Month)

ppt_lags2 <- ppt2 %>%
  arrange(Year, season) %>%
  #this creates a column for every lag this season to 5 ago
  do(data.frame(., setNames(shift(.$PPT, 1:5), c("PPT_l1", 'PPT_l2', "PPT_l3",
                                                  "PPT_l4", "PPT_l5")))) %>%
  ungroup() %>%
  #start in october, which is end of wet, 
  #so we can start with "wet"
  filter(season == 2) %>%
  dplyr::select(-season)

#get npp lags:
npp_lags <- npp3 %>%
  group_by(site) %>%
  arrange(year, seasonnum)  %>%
  #this creates a column for every lag this year to 5 seasons ago
  do(data.frame(., setNames(shift(.$NPP, 1:10), c("NPP_l1", 'NPP_l2', "NPP_l3",
                                                 "NPP_l4", "NPP_l5",
                                                 "NPP_l6", 'NPP_l7', "NPP_l8",
                                                 "NPP_l9", "NPP_l10")))) %>%
  ungroup() %>%
  #consider "this season" e.g. fall - and then previous seasons from
  #that for each web
  filter(seasonnum == 2) %>%
  dplyr::select(-seasonnum) %>%
  rename("ecosystem" = "site") %>%
  mutate(site = case_when(ecosystem == "core_black" ~ "BOER",
                          ecosystem == "core_creosote" ~ "LATR"))


# Combine datasets --------------------------------------------------------

all_data <- stability2 %>%
  left_join(temp_lags2, by = c("YEAR" = "Year")) %>%
  left_join(ppt_lags2, by = c("YEAR" = "Year")) %>%
  left_join(npp_lags, by = c("site", "web", "YEAR" = "year"))

# Prep observed data to bind ----------------------------------------------

raw_obs2 <- raw_dat %>%
  pivot_wider(names_from = "type",
              values_from = "bray")

all_data2 <- all_data %>%
  left_join(raw_obs2, by = c("yrID", "siteID"))

# Check for correlation ---------------------------------------------------

library(ggcorrplot)

dat1 <- all_data %>%
  dplyr::select(Temp:Temp_l5, PPT:PPT_l5)

ggcorrplot(cor(dat1, use = "complete.obs"), 
           type = "upper", lab= TRUE)

dat2 <- all_data %>%
  dplyr::select(Temp:Temp_l5, NPP:NPP_l5)

ggcorrplot(cor(dat2, use = "complete.obs"), 
           type = "upper", lab= TRUE)

dat3 <- all_data %>%
  dplyr::select(PPT:PPT_l5, NPP:NPP_l5)

ggcorrplot(cor(dat3, use = "complete.obs"), 
           type = "upper", lab = TRUE)
# Export ------------------------------------------------------------------

write.csv(all_data2,
          here("03_sev_grasshoppers",
               "data_outputs",
               'SAM',
               'data_prep',
               "sev_stability_metrics_with_covariates.csv"))





