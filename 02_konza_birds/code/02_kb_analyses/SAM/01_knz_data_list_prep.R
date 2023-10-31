#Prepping data object for the stability SAM models
#Ana Miller-ter Kuile
#July 5, 2023

#this is a script that preps data list for the turnover SAM models

# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 

package.list <- c("here", "tidyverse",
                  "data.table")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

## And loading them

for(i in package.list){library(i, character.only = T)}


# Load data ---------------------------------------------------------------

all_data <- read.csv(here('02_konza_birds',
                          'data_outputs',
                          "SAM",
                          "data_prep",
                          "stability_metrics_with_covariates.csv"))


# Remove years before 1984 ------------------------------------------------
all_data <- all_data %>%
  filter(RECYEAR >1983)
# Prep data for jags ------------------------------------------------------



n.data <- nrow(all_data)

bray <- as.vector(all_data$bray)
var.estimate <- as.vector(all_data$SD^2)

n.templag <- all_data %>%
  dplyr::select(TAVE:TAVE_l5) %>%
  ncol()

Temp <- all_data %>%
  dplyr::select(yrID, TransID, TAVE:TAVE_l5) %>%
  pivot_longer(TAVE:TAVE_l5,
               names_to = 'lag',
               values_to = 'temp') %>%
  mutate(temp = scale(temp)) %>%
  pivot_wider(names_from = 'lag',
              values_from = "temp") %>%
  dplyr::select(TAVE:TAVE_l5) %>%
  as.matrix()

sum(is.na(Temp))/(sum(is.na(Temp)) + sum(!is.na(Temp)))

PPT <- all_data %>%
  dplyr::select(yrID, TransID, PPT:PPT_l5) %>%
  pivot_longer(PPT:PPT_l5,
               names_to = 'lag',
               values_to = 'ppt') %>%
  mutate(ppt = scale(ppt)) %>%
  pivot_wider(names_from = 'lag',
              values_from = "ppt") %>%
  dplyr::select(PPT:PPT_l5) %>%
  as.matrix()

sum(is.na(PPT))/(sum(is.na(PPT)) + sum(!is.na(PPT)))

#< 4% missing data for PPT and TEmp

n.npplag <- all_data %>%
  dplyr::select(NPP:NPP_l5) %>%
  ncol()

NPP <- all_data %>%
  dplyr::select(yrID, TransID, NPP:NPP_l5) %>%
  pivot_longer(NPP:NPP_l5,
               names_to = 'lag',
               values_to = 'NPP') %>%
  mutate(NPP = scale(NPP)) %>%
  pivot_wider(names_from = 'lag',
              values_from = "NPP") %>%
  dplyr::select(NPP:NPP_l5) %>%
  as.matrix()

sum(is.na(NPP))/(sum(is.na(NPP)) + sum(!is.na(NPP)))

#16% missing data for NPP
#with subset data, 10% missing

# Make data list ----------------------------------------------------------

data <- list(n.data = n.data,
             bray = bray,
             var.estimate = var.estimate,
             n.templag = n.templag,
             Temp = Temp,
             PPT = PPT,
             n.npplag = n.npplag,
             NPP = NPP)

saveRDS(data, here('02_konza_birds',
                   "data_outputs",
                   'SAM',
                   "model_inputs",
                   "knz_bray_SAM_input_data.RDS"))


