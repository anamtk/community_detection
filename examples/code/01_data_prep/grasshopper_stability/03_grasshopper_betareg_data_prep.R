#Sevilleta grasshoppers environmental variable prep
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

#these are large files, so you can download them from here:
#https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-sev.1.16

c1 <- read.csv(here('examples',
                    'data_raw',
                    'grasshopper_stability',
                    'Sevilleta_LTER_Hourly_Meteorological_Data_1995_1999.csv'))

c2 <- read.csv(here('examples',
                    'data_raw',
                    'grasshopper_stability',
                    'Sevilleta_LTER_Hourly_Meteorological_Data_2000_2004.csv'))

c3 <- read.csv(here('examples',
                    'data_raw',
                    'grasshopper_stability',
                    'Sevilleta_LTER_Hourly_Meteorological_Data_2005_2009.csv'))

c4 <- read.csv(here('examples',
                    'data_raw',
                    'grasshopper_stability',
                    'Sevilleta_LTER_Hourly_Meteorological_Data_2010_2014.csv'))

c5 <- read.csv(here('examples',
                    'data_raw',
                    'grasshopper_stability',
                    'Sevilleta_LTER_Hourly_Meteorological_Data_2015_2019.csv'))

c6 <- read.csv(here('examples',
                    'data_raw',
                    'grasshopper_stability',
                    'Sevilleta_LTER_Hourly_Meteorological_Data_2020_2022.csv'))

#found here:
#https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-sev.331.2
npp <- read.csv(here('examples',
                     'data_raw',
                     'grasshopper_stability',
                     'sev331_quadrat_plant_species_biomass.csv'))

#to get the sites and transects we need for the other 
#two datasets
stability <- readRDS(here('examples',
                          'data_output',
                          'grasshopper_stability',
                          'tidy_data',
                          "grasshopper_stability_meanSD_corrected.RDS"))

IDs <- read.csv(here('examples',
                     'data_output',
                     'grasshopper_stability',
                     'tidy_data',
                     'grasshopper_msam_site_year_IDs.csv'))

raw_dat <- readRDS(here('examples',
                        'data_output',
                        'grasshopper_stability',
                        'tidy_data',
                        'grasshopper_stability_raw.RDS'))

# Get response data in order ----------------------------------------------

stability2 <- stability %>%
  dplyr::select(yrID, siteID, mean_bray, sd_bray,
                mean_balanced, sd_balanced, 
                mean_gradient, sd_gradient) %>%
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
  left_join(raw_obs2, by = c("yrID", "siteID")) %>%
  ungroup()

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
          here('examples',
               'data_output',
               'grasshopper_stability',
               'tidy_data',
               "grasshopper_betareg_tidydata.csv"))


# Prep data for jags ------------------------------------------------------


# prep modeled data for jags ----------------------------------------------

# Filter years with data on climate/npp -----------------------------------

all_data3 <- all_data2 %>%
  filter(YEAR > 1998) %>%
  #make site-web ID
  unite(c(site, web),
        col = "site_web",
        remove = F,
        sep = "_") %>%
  mutate(Web.ID = as.numeric(as.factor(site_web)))

# Prep data objects for model ---------------------------------------------


# Loop indexing -----------------------------------------------------------

n.data <- nrow(all_data3)

n.webs <- length(unique(all_data3$site_web))

n.transects <- length(unique(all_data3$site_web_trans))

n.templag <- all_data3 %>%
  dplyr::select(Temp:Temp_l5) %>%
  ncol()

n.pptlag <- all_data3 %>%
  dplyr::select(PPT:PPT_l5) %>%
  ncol()

n.npplag <- all_data3 %>%
  dplyr::select(NPP:NPP_l10) %>%
  ncol()

# Response data -----------------------------------------------------------

bray <- as.vector(all_data3$mean_bray)
var.estimate <- as.vector(all_data3$sd_bray^2)

# Random effects ----------------------------------------------------------

# Web.ID <- all_data3 %>%
#   distinct(siteID, Web.ID) %>%
#   dplyr::select(Web.ID) %>%
#   as_vector()

Web.ID <- as.vector(all_data3$Web.ID)

Transect.ID <- as.vector(all_data3$siteID)

# Covariates --------------------------------------------------------------

Temp <- all_data3  %>%
  dplyr::select(site_web_trans, YEAR, Temp:Temp_l5) %>%
  pivot_longer(Temp:Temp_l5,
               names_to = 'lag',
               values_to = 'temp') %>%
  mutate(temp = scale(temp)) %>%
  pivot_wider(names_from = 'lag',
              values_from = "temp") %>%
  dplyr::select(Temp:Temp_l5) %>%
  as.matrix()

sum(is.na(Temp))/(sum(is.na(Temp)) + sum(!is.na(Temp)))
#7% missing data
#with more lags - 0.4%

PPT <- all_data3  %>%
  dplyr::select(site_web_trans, YEAR, PPT:PPT_l5) %>%
  pivot_longer(PPT:PPT_l5,
               names_to = 'lag',
               values_to = 'ppt') %>%
  mutate(ppt = scale(ppt)) %>%
  pivot_wider(names_from = 'lag',
              values_from = "ppt") %>%
  dplyr::select(PPT:PPT_l5) %>%
  as.matrix()

sum(is.na(PPT))/(sum(is.na(PPT)) + sum(!is.na(PPT)))
#<5% missing data

NPP <- all_data3 %>%
  dplyr::select(site_web_trans, YEAR, NPP:NPP_l10) %>%
  pivot_longer(NPP:NPP_l10,
               names_to = 'lag',
               values_to = 'npp') %>%
  mutate(npp = scale(npp)) %>%
  pivot_wider(names_from = 'lag',
              values_from = "npp") %>%
  dplyr::select(NPP:NPP_l10) %>%
  as.matrix()

sum(is.na(NPP))/(sum(is.na(NPP)) + sum(!is.na(NPP)))
#~10% missing data


# Combine data into a data list -------------------------------------------


data1 <- list(n.data = n.data,
             n.webs = n.webs,
             n.templag = n.templag,
             n.npplag = n.npplag,
             n.pptlag = n.pptlag,
             n.transects = n.transects,
             bray = bray,
             var.estimate = var.estimate,
             Web.ID = Web.ID,
             Transect.ID = Transect.ID,
             Temp = Temp,
             PPT = PPT, 
             NPP = NPP)

saveRDS(data1, here('examples',
                   'data_output',
                   'grasshopper_stability',
                   'model_inputs',
                   'grasshopper_betareg_input_data_list_corrected.RDS'))



# Prep empirical data for jags --------------------------------------------

# Filter years with data on climate/npp -----------------------------------

all_data4 <- all_data2 %>%
  filter(YEAR > 1998) %>%
  #make site-web ID
  unite(c(site, web),
        col = "site_web",
        remove = F,
        sep = "_") %>%
  mutate(Web.ID = as.numeric(as.factor(site_web))) %>%
  filter(!is.na(observed_all))

# Prep data objects for model ---------------------------------------------


# Loop indexing -----------------------------------------------------------

n.data <- nrow(all_data4)

n.webs <- length(unique(all_data4$site_web))

n.transects <- length(unique(all_data4$site_web_trans))

n.templag <- all_data4 %>%
  dplyr::select(Temp:Temp_l5) %>%
  ncol()

n.pptlag <- all_data4 %>%
  dplyr::select(PPT:PPT_l5) %>%
  ncol()

n.npplag <- all_data4 %>%
  dplyr::select(NPP:NPP_l10) %>%
  ncol()

# Response data -----------------------------------------------------------

bray <- as.vector((all_data4$observed_all*(n.data-1) + 0.5)/n.data)


# Random effects ----------------------------------------------------------

# Web.ID <- all_data4 %>%
#   distinct(siteID, Web.ID) %>%
#   dplyr::select(Web.ID) %>%
#   as_vector()

Transect.ID <- as.vector(all_data4$siteID)

Web.ID <- as.vector(all_data4$Web.ID)

# Covariates --------------------------------------------------------------

Temp <- all_data4  %>%
  dplyr::select(site_web_trans, YEAR, Temp:Temp_l5) %>%
  pivot_longer(Temp:Temp_l5,
               names_to = 'lag',
               values_to = 'temp') %>%
  mutate(temp = scale(temp)) %>%
  pivot_wider(names_from = 'lag',
              values_from = "temp") %>%
  dplyr::select(Temp:Temp_l5) %>%
  as.matrix()

sum(is.na(Temp))/(sum(is.na(Temp)) + sum(!is.na(Temp)))
#7% missing data
#with more lags - 0.4%

PPT <- all_data4  %>%
  dplyr::select(site_web_trans, YEAR, PPT:PPT_l5) %>%
  pivot_longer(PPT:PPT_l5,
               names_to = 'lag',
               values_to = 'ppt') %>%
  mutate(ppt = scale(ppt)) %>%
  pivot_wider(names_from = 'lag',
              values_from = "ppt") %>%
  dplyr::select(PPT:PPT_l5) %>%
  as.matrix()

sum(is.na(PPT))/(sum(is.na(PPT)) + sum(!is.na(PPT)))
#<5% missing data

NPP <- all_data4 %>%
  dplyr::select(site_web_trans, YEAR, NPP:NPP_l10) %>%
  pivot_longer(NPP:NPP_l10,
               names_to = 'lag',
               values_to = 'npp') %>%
  mutate(npp = scale(npp)) %>%
  pivot_wider(names_from = 'lag',
              values_from = "npp") %>%
  dplyr::select(NPP:NPP_l10) %>%
  as.matrix()

sum(is.na(NPP))/(sum(is.na(NPP)) + sum(!is.na(NPP)))
#~10% missing data


# Combine data into a data list -------------------------------------------


data <- list(n.data = n.data,
             n.webs = n.webs,
             n.templag = n.templag,
             n.npplag = n.npplag,
             n.pptlag = n.pptlag,
             n.transects = n.transects,
             bray = bray,
             #var.estimate = var.estimate,
             Web.ID = Web.ID,
             Transect.ID = Transect.ID,
             Temp = Temp,
             PPT = PPT, 
             NPP = NPP)

saveRDS(data, here('examples',
                   'data_output',
                   'grasshopper_stability',
                   'model_inputs',
                   'grasshopper_betareg_input_data_list_raw.RDS'))











