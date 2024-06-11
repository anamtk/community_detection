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

all_data <- read.csv(here("03_sev_grasshoppers",
                          "data_outputs",
                          'SAM',
                          'data_prep',
                          "sev_stability_metrics_with_covariates.csv"))


# Filter years with data on climate/npp -----------------------------------

all_data2 <- all_data %>%
  filter(YEAR > 1998) %>%
  #make site-web ID
  unite(c(site, web),
        col = "site_web",
        remove = F,
        sep = "_") %>%
  mutate(Web.ID = as.numeric(as.factor(site_web)))

# Prep data objects for model ---------------------------------------------


# Loop indexing -----------------------------------------------------------

n.data <- nrow(all_data2)

n.webs <- length(unique(all_data2$site_web))

n.transects <- length(unique(all_data2$site_web_trans))

n.templag <- all_data2 %>%
  dplyr::select(Temp:Temp_l5) %>%
  ncol()

n.pptlag <- all_data2 %>%
  dplyr::select(PPT:PPT_l5) %>%
  ncol()

n.npplag <- all_data2 %>%
  dplyr::select(NPP:NPP_l10) %>%
  ncol()

# Response data -----------------------------------------------------------

bray <- as.vector(all_data2$bray)
var.estimate <- as.vector(all_data2$SD^2)

# Random effects ----------------------------------------------------------

# Web.ID <- all_data2 %>%
#   distinct(siteID, Web.ID) %>%
#   dplyr::select(Web.ID) %>%
#   as_vector()

Web.ID <- as.vector(all_data2$Web.ID)

Transect.ID <- as.vector(all_data2$siteID)

# Covariates --------------------------------------------------------------

Temp <- all_data2  %>%
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

PPT <- all_data2  %>%
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

NPP <- all_data2 %>%
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
             var.estimate = var.estimate,
             Web.ID = Web.ID,
             Transect.ID = Transect.ID,
             Temp = Temp,
             PPT = PPT, 
             NPP = NPP)

saveRDS(data, here('03_sev_grasshoppers',
                   "data_outputs",
                   'SAM',
                   "model_inputs",
                   "sev_bray_SAM_input_data_notrans.RDS"))









