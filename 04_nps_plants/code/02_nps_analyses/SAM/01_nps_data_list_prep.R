
#Prepping data object for the stability SAM models
#November 15, 2023

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

all_data <- read.csv(here("04_nps_plants",
                          "data_outputs",
                          'SAM',
                          'data_prep',
                          "nps_stability_metrics_with_covariates.csv"))


# Filter years with data on climate/npp -----------------------------------

all_data2 <- all_data %>%
  filter(!is.na(mean)) %>%
  #there should be 70 unique quadrat IDs, I believe
  unite(c(Plot, Transect, Quadrat),
        col = "quadrat_num",
        remove = F,
        sep = "_") %>%
#I think this is wrong:
  #mutate(Quad.ID = as.numeric(as.factor(quadnum))) %>%
  #make transect_num
  unite(c(Plot, Transect),
        col = "transect_num",
        remove = F,
        sep = "_") %>%
  mutate(Transect.ID = as.numeric(as.factor(transect_num))) %>%
  mutate(Plot.ID = as.numeric(as.factor(Plot))) %>%
  mutate(Quad.ID = as.numeric(as.factor(quadrat_num)))

# Prep data objects for model ---------------------------------------------


# Loop indexing -----------------------------------------------------------

n.data <- nrow(all_data2)

#I think there should be 70 unique quadrats
n.quads <- length(unique(all_data2$quadrat_num))

n.transects <- length(unique(all_data2$transect_num))

n.plots <- length(unique(all_data2$Plot))

n.lag <- all_data2 %>%
  dplyr::select(PPT:PPT_l20) %>%
  ncol()


# Response data -----------------------------------------------------------

# diss <- as.vector(all_data2$mean)
# var.estimate <- as.vector(all_data2$sd)

diss <- as.vector(all_data2$mean_gain)
var.estimate <- as.vector(all_data2$sd_gain)

# Random effects ----------------------------------------------------------

Plot.ID  <- all_data2 %>%
  distinct(Transect.ID, Plot.ID) %>%
  dplyr::select(Plot.ID) %>%
  as_vector()

Transect.ID <- as.vector(all_data2$Transect.ID)

Quad.ID <- as.vector(all_data2$Quad.ID)

# Covariates --------------------------------------------------------------

PPT <- all_data2  %>%
  dplyr::select(quadnum, EventYear, PPT:PPT_l20) %>%
  pivot_longer(PPT:PPT_l20,
               names_to = 'lag',
               values_to = 'ppt') %>%
  mutate(ppt = scale(ppt)) %>%
  pivot_wider(names_from = 'lag',
              values_from = "ppt") %>%
  dplyr::select(PPT:PPT_l20) %>%
  as.matrix()

sum(is.na(PPT))/(sum(is.na(PPT)) + sum(!is.na(PPT)))
#~8% missing data

VPD <- all_data2 %>%
  dplyr::select(quadnum, EventYear, VPD:VPD_l20) %>%
  pivot_longer(VPD:VPD_l20,
               names_to = 'lag',
               values_to = 'vpd') %>%
  mutate(vpd = scale(vpd)) %>%
  pivot_wider(names_from = 'lag',
              values_from = "vpd") %>%
  dplyr::select(VPD:VPD_l20) %>%
  as.matrix()

sum(is.na(VPD))/(sum(is.na(VPD)) + sum(!is.na(VPD)))
#~8% missing data


# Combine data into a data list -------------------------------------------


data <- list(n.data = n.data,
             n.lag = n.lag,
             n.quads = n.quads,
             n.transects = n.transects,
             n.plots = n.plots,
             diss = diss,
             var.estimate = var.estimate,
             Plot.ID = Plot.ID,
             Transect.ID = Transect.ID,
             Quad.ID = Quad.ID,
             PPT = PPT, 
             VPD = VPD)

saveRDS(data, here('04_nps_plants',
                   "data_outputs",
                   'SAM',
                   "model_inputs",
                   "nps_diss_SAM_input_data.RDS"))

