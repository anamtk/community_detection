
#Prepping data object for the stability SAM models
#November 15, 2023

#this is a script that preps data list for the turnover SAM models

# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 

package.list <- c("here", "tidyverse",
                  "data.table","jagsUI",
                  'rjags',
                  'mcmcplots',
                  "coda",
                  'patchwork')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

## And loading them

for(i in package.list){library(i, character.only = T)}


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


# Gains -------------------------------------------------------------------

gaindf <- all_data2 %>%
  filter(!is.na(mean_gain))

# Prep data objects for model ---------------------------------------------


# Loop indexing -----------------------------------------------------------

n.data <- nrow(gaindf)

#I think there should be 70 unique quadrats
n.quads <- length(unique(gaindf$quadrat_num))

n.transects <- length(unique(gaindf$transect_num))

n.plots <- length(unique(gaindf$Plot))

n.lag <- gaindf %>%
  dplyr::select(PPT:PPT_l7) %>%
  ncol()


# Response data -----------------------------------------------------------

# diss <- as.vector(gaindf2$mean)
# var.estimate <- as.vector(gaindf2$sd^2)

diss <- as.vector(gaindf$mean_gain)
var.estimate <- as.vector(gaindf$sd_gain^2)


# Random effects ----------------------------------------------------------

Plot.ID  <- gaindf %>%
  distinct(Transect.ID, Plot.ID) %>%
  dplyr::select(Plot.ID) %>%
  as_vector()

Transect.ID <- as.vector(gaindf$Transect.ID)

Quad.ID <- as.vector(gaindf$Quad.ID)

# Covariates --------------------------------------------------------------

PPT <- gaindf  %>%
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
#~0.6% missing data (new - 5/6/24)
#~8% missing data (old)

VPD <- gaindf %>%
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
#~0.6% missing data (new - 5/6/24)
#~8% missing data (old)


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


# Run model ---------------------------------------------------------------

# Parameters to save ------------------------------------------------------

params <- c(#'b0.quad',
  'b0.transect',
  #'b0.plot',
  'b',
  'b0',
  'wA',
  'wB',
  #'sig.quad',
  'sig.transect',
  #'sig.plot',
  'var.process')



# JAGS model --------------------------------------------------------------

model <- here('examples',
              'code',
              '02_models',
              '02_beta_regression',
              'plant_turnover',
              'detection_corrected',
              'plant_betareg_model_corrected.R')

Sys.time()
mod <- jagsUI::jags(data = data,
                    inits = NULL,
                    model.file = model,
                    parameters.to.save = params,
                    parallel = TRUE,
                    n.chains = 3,
                    n.burnin = 1000,
                    n.iter =  40000,
                    n.thin = 7,
                    DIC = TRUE)

Sys.time()

# # Check convergence -------------------------------------------------------
# 
mcmcplot(mod$samples)
# 
gelman.diag(mod$samples, multivariate = F)
# 
#
sum <- summary(mod$samples)

saveRDS(sum, here('examples',
                  'data_output',
                  'plant_turnover',
                  'model_outputs',
                  'plant_gains_betareg_summary_corrected.RDS'))


# Losses ------------------------------------------------------------------

lossdf <- all_data2 %>%
  filter(!is.na(mean_loss))

# Prep data objects for model ---------------------------------------------


# Loop indexing -----------------------------------------------------------

n.data <- nrow(lossdf)

#I think there should be 70 unique quadrats
n.quads <- length(unique(lossdf$quadrat_num))

n.transects <- length(unique(lossdf$transect_num))

n.plots <- length(unique(lossdf$Plot))

n.lag <- lossdf %>%
  dplyr::select(PPT:PPT_l7) %>%
  ncol()


# Response data -----------------------------------------------------------

# diss <- as.vector(lossdf$mean)
# var.estimate <- as.vector(lossdf$sd^2)

diss <- as.vector(lossdf$mean_loss)
var.estimate <- as.vector(lossdf$sd_loss^2)


# Random effects ----------------------------------------------------------

Plot.ID  <- lossdf %>%
  distinct(Transect.ID, Plot.ID) %>%
  dplyr::select(Plot.ID) %>%
  as_vector()

Transect.ID <- as.vector(lossdf$Transect.ID)

Quad.ID <- as.vector(lossdf$Quad.ID)

# Covariates --------------------------------------------------------------

PPT <- lossdf  %>%
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
#~0.6% missing data (new - 5/6/24)
#~8% missing data (old)

VPD <- lossdf %>%
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
#~0.6% missing data (new - 5/6/24)
#~8% missing data (old)


# Combine data into a data list -------------------------------------------


data2 <- list(n.data = n.data,
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


# Run model ---------------------------------------------------------------

# JAGS model --------------------------------------------------------------

Sys.time()
mod2 <- jagsUI::jags(data = data2,
                     inits = NULL,
                     model.file = model,
                     parameters.to.save = params,
                     parallel = TRUE,
                     n.chains = 3,
                     n.burnin = 1000,
                     n.iter =  40000,
                     n.thin = 7,
                     DIC = TRUE)

Sys.time()

# # Check convergence -------------------------------------------------------
# 
mcmcplot(mod2$samples)
# 
gelman.diag(mod2$samples, multivariate = F)
# 
#
sum2 <- summary(mod2$samples)

saveRDS(sum2, here('examples',
                   'data_output',
                   'plant_turnover',
                   'model_outputs',
                   'plant_loss_betareg_summary_corrected.RDS'))


