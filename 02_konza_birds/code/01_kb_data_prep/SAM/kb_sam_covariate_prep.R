# Covariate prep for SAM for Konza bird dataset
# Ana Miller-ter Kuile
# October 11, 2023

#this script preps covariates for the SAM model for Konza birds

# Load packages -----------------------------------------------------------


package.list <- c("here", "tidyverse",
                  'data.table')

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


# Get those monthly values to be in lag format for model ------------------

#get temp lags:
temp_lags <- temp %>%
  arrange(RECYEAR, RECMONTH) %>%
  #this creates a column for every lag this month to 11 months ago
  do(data.frame(., setNames(shift(.$TAVE, 1:11), c("TAVE_l1", 'TAVE_l2', "TAVE_l3",
                                                   "TAVE_l4", "TAVE_l5", "TAVE_l6",
                                                   "TAVE_l7", "TAVE_l8", "TAVE_l9",
                                                   "TAVE_l10", "TAVE_l11")))) %>%
  ungroup()

#get ppt lags:
ppt_lags <- ppt %>%
  arrange(RECYEAR, RECMONTH) %>%
  #this creates a column for every lag this month to 11 months ago
  do(data.frame(., setNames(shift(.$PPT, 1:11), c("PPT_l1", 'PPT_l2', "PPT_l3",
                                                   "PPT_l4", "PPT_l5", "PPT_l6",
                                                   "PPT_l7", "PPT_l8", "PPT_l9",
                                                   "PPT_l10", "PPT_l11")))) %>%
  ungroup()



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
  ungroup()


