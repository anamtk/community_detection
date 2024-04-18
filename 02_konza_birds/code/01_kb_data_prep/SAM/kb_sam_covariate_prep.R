# Covariate prep for SAM for Konza bird dataset
# Ana Miller-ter Kuile
# October 11, 2023

#this script preps covariates for the SAM model for Konza birds

# Load packages -----------------------------------------------------------


package.list <- c("here", "tidyverse",
                  'data.table', 'ggcorrplot')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load data ---------------------------------------------------------------

#we'll be using some summary of the weather station precip and temperatures
#(maybe monthly??) 
#and also plant biomass data from the plots - this is yearly I think, but there
#are lots of quadrats - so we can get an average
#daily climate
climate <- read.csv(here('02_konza_birds',
                         'data_raw',
                         'environmental',
                         'AWE012.csv'))

#yearly npp for a bunch of quadrats in each watershed
npp <- read.csv(here('02_konza_birds',
                     'data_raw',
                     'environmental',
                     'PAB011.csv'))


#to get the sites and transects we need for the other 
  #two datasets
data <- readRDS(here("02_konza_birds",
                     "monsoon",
                     "MSAM",
                     "outputs",
                     "bird_bray_meanSD.RDS"))

IDs <- read.csv(here('02_konza_birds',
                     "data_outputs",
                     "metadata",
                     "site_year_IDs.csv"))

konza_obs <- readRDS(here('05_visualizations',
                          'viz_data',
                          'konza_observed_bray.RDS'))

# Summarize monthly climate values ----------------------------------------

#average mean temp for the month
#average precipitation for the month

temp <- climate %>%
  dplyr::select(RECYEAR, RECMONTH, TAVE) %>%
  mutate(TAVE = as.numeric(TAVE)) %>%
  group_by(RECYEAR, RECMONTH) %>%
  summarise(TAVE = mean(TAVE, na.rm = T)) %>%
  ungroup()

ppt <- climate %>%
  dplyr::select(RECYEAR, RECMONTH, DPPT) %>%
  mutate(DPPT = case_when(DPPT == "." ~"0",
                          TRUE ~ DPPT)) %>%
  mutate(DPPT = as.numeric(DPPT)) %>%
  group_by(RECYEAR, RECMONTH) %>%
  summarise(PPT = mean(DPPT, na.rm = T)) %>%
  ungroup()

ggplot(temp, aes(x = RECMONTH, y = TAVE)) +
  geom_point(position = position_jitter(width = 0.25))

ggplot(ppt, aes(x = RECMONTH, y = PPT)) +
  geom_point(position = position_jitter(width = 0.25))

#standardize:
#whatever months have + are one season, whatever months
#have - values are the other season
#1,2,3
#11, 12
temp %>%
  mutate(TAVE = scale(TAVE)) %>%
  group_by(RECMONTH) %>%
  summarise(mean = mean(TAVE))

#1, 2, 3
#10, 11, 12
ppt %>%
  mutate(PPT = scale(PPT)) %>%
  group_by(RECMONTH) %>%
  summarise(mean = mean(PPT))

#going to summarise both by the same months - using the PPT values
#since these split things nicely 6-6, and month 10 is *barely* above
#0 in the temp dataset

temp2 <- climate %>%
  dplyr::select(RECYEAR, RECMONTH, TAVE) %>%
  mutate(TAVE = as.numeric(TAVE)) %>%
  mutate(season = case_when(RECMONTH %in% c(1, 2, 3, 10, 11, 12) ~ "winter",
                            RECMONTH %in% c(4, 5, 6, 7, 8, 9) ~ "summer",
                            TRUE ~ NA_character_)) %>%
  group_by(RECYEAR, season) %>%
  summarise(TAVE = mean(TAVE, na.rm = T)) %>%
  ungroup()

ppt2 <- climate %>%
  dplyr::select(RECYEAR, RECMONTH, DPPT) %>%
  mutate(DPPT = as.numeric(DPPT)) %>%
  mutate(season = case_when(RECMONTH %in% c(1, 2, 3, 10, 11, 12) ~ "winter",
                            RECMONTH %in% c(4, 5, 6, 7, 8, 9) ~ "summer",
                            TRUE ~ NA_character_)) %>%
  group_by(RECYEAR, season) %>%
  summarise(PPT = mean(DPPT, na.rm = T)) %>%
  ungroup()


# Get seasonal values into lag format for model ---------------------------

#since sampling usually occured in June and January, but sometimes
#in February, April, October, and December,
#going to consider lags to be "end of year and back",
#so starting in "winter" of each year and going back in time
#get temp lags:
temp_lags <- temp2 %>%
  arrange(RECYEAR, season) %>%
  #this creates a column for every lag this season (winter) and 5 seasons ago
  do(data.frame(., setNames(shift(.$TAVE, 1:5), c("TAVE_l1", 'TAVE_l2', "TAVE_l3",
                                                   "TAVE_l4", "TAVE_l5")))) %>%
  ungroup() %>%
  filter(season == "winter") %>%
  dplyr::select(-season)

#get ppt lags:
ppt_lags <- ppt2 %>%
  arrange(RECYEAR, season) %>%
  #this creates a column for every lag this season (winter) and 5 seasons ago
  do(data.frame(., setNames(shift(.$PPT, 1:11), c("PPT_l1", 'PPT_l2', "PPT_l3",
                                                   "PPT_l4", "PPT_l5", "PPT_l6",
                                                  "PPT_l7", "PPT_l8", "PPT_l9",
                                                  "PPT_l10", "PPT_l11")))) %>%
  ungroup() %>%
  filter(season == "winter") %>%
  dplyr::select(-season)



# Summarise NPP data ------------------------------------------------------

npp2 <- npp %>%
  rowwise() %>%
  mutate(totallive = sum(LVGRASS, FORBS, na.rm = T)) %>%
  group_by(RECYEAR, WATERSHED) %>%
  summarise(NPP = mean(totallive, na.rm = T)) %>%
  ungroup()


# Set up NPP to be lags ---------------------------------------------------

#get npp lags:
npp_lags <- npp2 %>%
  group_by(WATERSHED) %>%
  arrange(RECYEAR) %>%
  #this creates a column for every lag this year to 5 years ago
  do(data.frame(., setNames(shift(.$NPP, 1:5), c("NPP_l1", 'NPP_l2', "NPP_l3",
                                                  "NPP_l4", "NPP_l5")))) %>%
  ungroup() %>%
  mutate(WATERSHED = case_when(WATERSHED == "001d" ~ "001D",
                               WATERSHED == "004b" ~ "004B",
                               WATERSHED == "020b" ~ "020B",
                               TRUE ~ NA_character_))


# Prep stability data -----------------------------------------------------

stability2 <- as.data.frame(data) %>%
  rownames_to_column(var = "var") %>%
  filter(var != "deviance") %>%
  separate(var,
           into = c('TransID', 'yrID'),
           sep = ",") %>%
  mutate(TransID = str_sub(TransID, 6, nchar(TransID))) %>%
  mutate(yrID = str_sub(yrID, 1, (nchar(yrID)-1))) %>%
  rename("bray" = "Mean") %>%
  dplyr::select(yrID, TransID, bray, SD) %>%
  mutate(yrID = as.numeric(yrID),
         TransID = as.numeric(TransID)) %>%
  left_join(IDs, by = c("TransID", "yrID")) %>%
  dplyr::select(-X)

# Get site and transect IDs -----------------------------------------------

sites <- stability2 %>%
  distinct(TransID, TRANSNUM) 

site <- sites$SITE
sitetrans <- sites$SITE_TRANS

# Combine all data --------------------------------------------------------

all_data <- stability2 %>%
  left_join(npp_lags, by = c("RECYEAR",
                             "WATERSHED")) %>%
  left_join(temp_lags, by = c("RECYEAR")) %>%
  left_join(ppt_lags, by = c("RECYEAR")) 


# Prep observed data to bind ----------------------------------------------

konza_obs2 <- konza_obs %>%
  pivot_wider(names_from = "type",
              values_from = "bray")

all_data2 <- all_data %>%
  left_join(konza_obs2, by = c("yrID", "TransID"))

# Check for correlation in predictors -------------------------------------

corr_check <- all_data %>%
  dplyr::select(NPP:PPT_l5) %>%
  filter(if_all(everything(), ~!is.na(.))) 
  
corr_check2 <- cor(corr_check)
ggcorrplot(cor(corr_check), type = "lower", lab = T)

# Export ------------------------------------------------------------------



write.csv(all_data2, here("02_konza_birds",
                         "data_outputs",
                         'SAM',
                         'data_prep',
                         "stability_metrics_with_covariates.csv"))

