
#Prepping data object for the stability SAM models
#November 15, 2023

#this is a script that preps data list for the turnover SAM models

# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 

package.list <- c("here", "tidyverse",
                  "data.table",
                  "jagsUI",
                  'rjags',
                  'mcmcplots',
                  "coda",
                  'patchwork')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

## And loading them

for(i in package.list){library(i, character.only = T)}


source(here('examples',
            'code',
            '00_functions',
            'plot_functions.R'))
# Load data ---------------------------------------------------------------

all_data <- read.csv(here('examples',
                          'data_output',
                          'plant_turnover',
                          'tidy_data',
                          'plant_betareg_tidydata.csv'))
                          
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


# Gain --------------------------------------------------------------------

gains <- all_data2 %>%
  filter(!is.na(gain))

# Loop indexing -----------------------------------------------------------

n.data <- nrow(gains)

#I think there should be 70 unique quadrats
n.quads <- length(unique(gains$quadrat_num))

n.transects <- length(unique(gains$transect_num))

n.plots <- length(unique(gains$Plot))

#going back 2 years
n.lag <- gains %>%
  dplyr::select(PPT:PPT_l7) %>%
  ncol()

# n.lag <- all_data2 %>%
#   dplyr::select(PPT:PPT_l20) %>%
#   ncol()


# Response data -----------------------------------------------------------

# diss <- as.vector(all_data2$mean)
# var.estimate <- as.vector(all_data2$sd^2)
diss <- as.vector(gains$gain)
diss <- as.vector((gains$gain*(n.data-1) + 0.5)/n.data)

# Random effects ----------------------------------------------------------

Plot.ID  <- gains %>%
  distinct(Transect.ID, Plot.ID) %>%
  dplyr::select(Plot.ID) %>%
  as_vector()

Transect.ID <- as.vector(gains$Transect.ID)

# Transect.ID  <- all_data2 %>%
#   distinct(Quad.ID, Transect.ID) %>%
#   dplyr::select(Transect.ID) %>%
#   as_vector()

Quad.ID <- as.vector(gains$Quad.ID)

# Covariates --------------------------------------------------------------

PPT <- gains  %>%
  dplyr::select(quadnum, EventYear, PPT:PPT_l7) %>%
  pivot_longer(PPT:PPT_l7,
               names_to = 'lag',
               values_to = 'ppt') %>%
  mutate(ppt = scale(ppt)) %>%
  pivot_wider(names_from = 'lag',
              values_from = "ppt") %>%
  dplyr::select(PPT:PPT_l7) %>%
  as.matrix()

sum(is.na(PPT))/(sum(is.na(PPT)) + sum(!is.na(PPT)))
#~8% missing data

VPD <- gains %>%
  dplyr::select(quadnum, EventYear, VPD:VPD_l7) %>%
  pivot_longer(VPD:VPD_l7,
               names_to = 'lag',
               values_to = 'vpd') %>%
  mutate(vpd = scale(vpd)) %>%
  pivot_wider(names_from = 'lag',
              values_from = "vpd") %>%
  dplyr::select(VPD:VPD_l7) %>%
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
             #var.estimate = var.estimate,
             Plot.ID = Plot.ID,
             Transect.ID = Transect.ID,
             Quad.ID = Quad.ID,
             PPT = PPT, 
             VPD = VPD)


# Run model ---------------------------------------------------------------

# Parameters to save ------------------------------------------------------

params <- c('b0.quad',
            'b0.transect',
            #'b0.plot',
            'b0',
            'b',
            'wA',
            'wB',
            'sig.quad',
            'sig.transect',
            #'sig.plot',
            'var.process')



# JAGS model --------------------------------------------------------------

model <- here('examples',
              'code',
              '02_models',
              '02_beta_regression',
              'plant_turnover',
              'uncorrected',
              'plant_betareg_model_raw.R')

Sys.time()
mod <- jagsUI::jags(data = data,
                    inits = NULL,
                    model.file = model,
                    parameters.to.save = params,
                    parallel = TRUE,
                    n.chains = 3,
                    n.burnin = 15000,
                    n.thin = 7, 
                    n.iter = 41000, #4000,
                    DIC = TRUE)

Sys.time()

# # Check convergence -------------------------------------------------------
# 
mcmcplot(mod$samples)
# 
gelman.diag(mod$samples, multivariate = F)
# 
rhat_graph_fun(list = mod$Rhat, parms = params, rhat = 1.1) +
  labs(title = "NPS plant SAM Rhat: \n observed data")



sum <- summary(mod$samples)

saveRDS(sum, here('examples',
                  'data_output',
                  'plant_turnover',
                  'model_outputs',
                  'plant_gains_betareg_summary_raw.RDS'))


# Losses ------------------------------------------------------------------

losses <- all_data2 %>%
  filter(!is.na(loss))

# Loop indexing -----------------------------------------------------------

n.data <- nrow(losses)

#I think there should be 70 unique quadrats
n.quads <- length(unique(losses$quadrat_num))

n.transects <- length(unique(losses$transect_num))

n.plots <- length(unique(losses$Plot))

#going back 2 years
n.lag <- losses %>%
  dplyr::select(PPT:PPT_l7) %>%
  ncol()

# n.lag <- all_data2 %>%
#   dplyr::select(PPT:PPT_l20) %>%
#   ncol()


# Response data -----------------------------------------------------------

# diss <- as.vector(all_data2$mean)
# var.estimate <- as.vector(all_data2$sd^2)
diss <- as.vector(losses$loss)
diss <- as.vector((losses$loss*(n.data-1) + 0.5)/n.data)

# Random effects ----------------------------------------------------------

Plot.ID  <- losses %>%
  distinct(Transect.ID, Plot.ID) %>%
  dplyr::select(Plot.ID) %>%
  as_vector()

Transect.ID <- as.vector(losses$Transect.ID)

# Transect.ID  <- all_data2 %>%
#   distinct(Quad.ID, Transect.ID) %>%
#   dplyr::select(Transect.ID) %>%
#   as_vector()

Quad.ID <- as.vector(losses$Quad.ID)

# Covariates --------------------------------------------------------------

PPT <- losses  %>%
  dplyr::select(quadnum, EventYear, PPT:PPT_l7) %>%
  pivot_longer(PPT:PPT_l7,
               names_to = 'lag',
               values_to = 'ppt') %>%
  mutate(ppt = scale(ppt)) %>%
  pivot_wider(names_from = 'lag',
              values_from = "ppt") %>%
  dplyr::select(PPT:PPT_l7) %>%
  as.matrix()

sum(is.na(PPT))/(sum(is.na(PPT)) + sum(!is.na(PPT)))
#~8% missing data

VPD <- losses %>%
  dplyr::select(quadnum, EventYear, VPD:VPD_l7) %>%
  pivot_longer(VPD:VPD_l7,
               names_to = 'lag',
               values_to = 'vpd') %>%
  mutate(vpd = scale(vpd)) %>%
  pivot_wider(names_from = 'lag',
              values_from = "vpd") %>%
  dplyr::select(VPD:VPD_l7) %>%
  as.matrix()

sum(is.na(VPD))/(sum(is.na(VPD)) + sum(!is.na(VPD)))
#~8% missing data


# Combine data into a data list -------------------------------------------


data2 <- list(n.data = n.data,
             n.lag = n.lag,
             n.quads = n.quads,
             n.transects = n.transects,
             n.plots = n.plots,
             diss = diss,
             #var.estimate = var.estimate,
             Plot.ID = Plot.ID,
             Transect.ID = Transect.ID,
             Quad.ID = Quad.ID,
             PPT = PPT, 
             VPD = VPD)

# saveRDS(data, here('04_nps_plants',
#                    "data_outputs",
#                    'SAM',
#                    "model_inputs",
#                    "nps_losses_SAM_input_data_raw.RDS"))
# JAGS model --------------------------------------------------------------

Sys.time()
mod2 <- jagsUI::jags(data = data2,
                     inits = NULL,
                     model.file = model,
                     parameters.to.save = params,
                     parallel = TRUE,
                     n.chains = 3,
                     n.burnin = 15000,
                     n.thin = 15,
                     n.iter = 50000, #4000,
                     DIC = TRUE)

Sys.time()

# # Check convergence -------------------------------------------------------
# 
mcmcplot(mod2$samples)
# 
gelman.diag(mod2$samples, multivariate = F)
# 
rhat_graph_fun(list = mod2$Rhat, parms = params, rhat = 1.1) +
  labs(title = "NPS plant SAM Rhat: \n observed data")

sum <- summary(mod$samples)

saveRDS(sum, here('examples',
                  'data_output',
                  'plant_turnover',
                  'model_outputs',
                  'plant_loss_betareg_summary_raw.RDS'))



