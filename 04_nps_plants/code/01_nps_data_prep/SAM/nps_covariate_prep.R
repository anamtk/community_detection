
#NPS plants environmental variable prep
# November 15, 2023

#this script preps the environmental variables with lags
#for the NPS Plants dataset


# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse", 'data.table')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Load data ---------------------------------------------------------------

climate <- read.csv(here('04_nps_plants',
                        'data_raw',
                        'seasonal_ppt_vpd.csv'))
climate <- climate[,-c(1:2)]

#to get the sites and transects we need for the other 
#two datasets
stability <- readRDS(here("04_nps_plants",
                          "monsoon",
                          "nps_MSAM",
                          "outputs_yrsite",
                          "nps_Jaccard_summary.RDS"))

IDs <- read.csv(here('04_nps_plants',
                     'data_outputs',
                     'metadata',
                     'site_year_IDs.csv'), row.names=1)


# Get response data in order ----------------------------------------------

stability2 <- as.data.frame(stability) %>%
  separate(plot_trans_quad,
           into = c('Plot', 'Transect','Quadrat'),
           sep = "_") %>%
  mutate(Quadrat = as.numeric(Quadrat)) %>%
  left_join(IDs, by = c("Plot", "Transect", "Quadrat", "EventYear"))




# ppt --------------------------------------------------------------------

ppt <- climate %>%
  dplyr::select(-c(vpd_monsoon:vpd_winter)) %>%
  pivot_longer(cols = ppt_monsoon:ppt_winter,
               names_to = "season",
               values_to = "PPT") %>%
  mutate(season = case_when(season == "ppt_spring" ~ 1,
                            season == "ppt_earlysummer" ~ 2,
                            season == "ppt_monsoon"~ 3,
                            TRUE ~ 4))

# vpd --------------------------------------------------------------------


vpd <- climate %>%
  dplyr::select(-c(ppt_monsoon:ppt_winter)) %>%
  pivot_longer(cols = vpd_monsoon:vpd_winter,
               names_to = "season",
               values_to = "VPD") %>%
  mutate(season = case_when(season == "vpd_spring" ~ 1,
                            season == "vpd_earlysummer" ~ 2,
                            season == "vpd_monsoon" ~ 3,
                            TRUE ~ 4))


# Get lags ----------------------------------------------------------------

ppt_lags <- ppt %>%
  group_by(Plot) %>%
  arrange(EventYear, season) %>%
  #this creates a column for every lag this season to 5 years ago
  do(data.frame(., setNames(shift(.$PPT, 1:20), c('PPT_l1', 'PPT_l2', 'PPT_l3',
                                                  'PPT_l4', 'PPT_l5',
                                                  'PPT_l6', 'PPT_l7',
                                                  'PPT_l8', 'PPT_l9',
                                                  'PPT_l10', 'PPT_l11', 
                                                  'PPT_l12', 'PPT_l13',
                                                  'PPT_l14','PPT_l15',
                                                  'PPT_l16','PPT_l17',
                                                  'PPT_l18','PPT_l19',
                                                  'PPT_l20')))) %>%
  ungroup() %>%
  #start in monsoon, which is when sampling ended
  filter(season == 3) %>%
  dplyr::select(-season)


#get npp lags:
vpd_lags <- vpd %>%
  group_by(Plot) %>%
  arrange(EventYear, season)  %>%
  #this creates a column for every lag this season to 5 years ago
  do(data.frame(., setNames(shift(.$VPD, 1:20), c('VPD_l1', 'VPD_l2', 'VPD_l3',
                                                  'VPD_l4', 'VPD_l5',
                                                  'VPD_l6', 'VPD_l7',
                                                  'VPD_l8', 'VPD_l9',
                                                  'VPD_l10', 'VPD_l11', 
                                                  'VPD_l12', 'VPD_l13',
                                                  'VPD_l14','VPD_l15',
                                                  'VPD_l16','VPD_l17',
                                                  'VPD_l18','VPD_l19',
                                                  'VPD_l20')))) %>%
  ungroup() %>%
  #start in monsoon, which is when sampling ended
  filter(season == 3) %>%
  dplyr::select(-season)


# Combine datasets --------------------------------------------------------

all_data <- stability2 %>%
  left_join(ppt_lags, by = c("Plot", "EventYear")) %>%
  left_join(vpd_lags, by = c("Plot", "EventYear"))


# Check for correlation ---------------------------------------------------

library(ggcorrplot)

dat1 <- all_data %>%
  dplyr::select(PPT:PPT_l20, VPD:VPD_l20)

ggcorrplot(cor(dat1, use = "complete.obs"), 
           type = "upper", lab= TRUE)


# Export ------------------------------------------------------------------

write.csv(all_data,
          here("04_nps_plants",
               "data_outputs",
               'SAM',
               'data_prep',
               "nps_stability_metrics_with_covariates.csv"))

