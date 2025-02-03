#Prepping data object for the stability SAM models
#Ana Miller-ter Kuile
#July 5, 2023

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

source(here('examples',
            'code',
            "00_functions",
            'plot_functions.R'))
# Load data ---------------------------------------------------------------

all_data <- read.csv(here('examples',
                          'data_output',
                          'bird_fundiv',
                          'tidy_data',
                          'bird_betareg_tidy_data.csv'))


# Remove years before 1982 ------------------------------------------------
all_data <- all_data %>%
  filter(RECYEAR >1981)

FRic <- all_data %>%
  filter(!is.na(FRic))

Q <- all_data %>%
  filter(!is.na(Q))

# Prep data for jags ------------------------------------------------------

# FRic --------------------------------------------------------------------

# Loop indexing -----------------------------------------------------------

n.data <- nrow(FRic)

n.templag <- FRic %>%
  dplyr::select(TAVE:TAVE_l3) %>%
  ncol()

n.npplag <- FRic %>%
  dplyr::select(NPP:NPP_l3) %>%
  ncol()

n.transects <- length(unique(FRic$TransID))


# Random effects ----------------------------------------------------------

Transect.ID <- FRic %>%
  dplyr::select(TransID) %>%
  as_vector()

# Response data -----------------------------------------------------------

#transform like: https://stats.stackexchange.com/questions/31300/dealing-with-0-1-values-in-a-beta-regression
#bray <- as.vector((all_data$observed_all*(n.data-1) + 0.5)/n.data)
fd <- as.vector((FRic$FRic))
#fd <- as.vector((FRic$FRic*(n.data-1) + 0.5)/n.data)

# Covariates --------------------------------------------------------------

Temp <- FRic %>%
  dplyr::select(yrID, TransID, TAVE:TAVE_l3) %>%
  pivot_longer(TAVE:TAVE_l3,
               names_to = 'lag',
               values_to = 'temp') %>%
  mutate(temp = scale(temp)) %>%
  pivot_wider(names_from = 'lag',
              values_from = "temp") %>%
  dplyr::select(TAVE:TAVE_l3) %>%
  as.matrix()

sum(is.na(Temp))/(sum(is.na(Temp)) + sum(!is.na(Temp)))

PPT <- FRic %>%
  dplyr::select(yrID, TransID, PPT:PPT_l3) %>%
  pivot_longer(PPT:PPT_l3,
               names_to = 'lag',
               values_to = 'ppt') %>%
  mutate(ppt = scale(ppt)) %>%
  pivot_wider(names_from = 'lag',
              values_from = "ppt") %>%
  dplyr::select(PPT:PPT_l3) %>%
  as.matrix()

sum(is.na(PPT))/(sum(is.na(PPT)) + sum(!is.na(PPT)))

# 0% missing data for PPT and TEmp

NPP <- FRic %>%
  dplyr::select(yrID, TransID, NPP:NPP_l3) %>%
  pivot_longer(NPP:NPP_l3,
               names_to = 'lag',
               values_to = 'NPP') %>%
  mutate(NPP = scale(NPP)) %>%
  pivot_wider(names_from = 'lag',
              values_from = "NPP") %>%
  dplyr::select(NPP:NPP_l3) %>%
  as.matrix()

sum(is.na(NPP))/(sum(is.na(NPP)) + sum(!is.na(NPP)))

#with subset data, 10% missing

# Make data list ----------------------------------------------------------

data <- list(n.data = n.data,
             n.templag = n.templag,
             n.npplag = n.npplag,
             n.transects = n.transects,
             fd = fd,
             Transect.ID = Transect.ID,
             Temp = Temp,
             PPT = PPT,
             NPP = NPP)

# Parameters to save ------------------------------------------------------

params <- c('b',
            'b0',
            'wA',
            'wB',
            #'wC',
            'var.process')



# JAGS model --------------------------------------------------------------

model <- here('examples',
              'code',
              '02_models',
              '02_beta_regression',
              'bird_fundiv',
              'uncorrected',
              'bird_betareg_model_raw.R')

Sys.time()
mod <- jagsUI::jags(data = data,
                    inits = NULL,
                    model.file = model,
                    parameters.to.save = params,
                    parallel = TRUE,
                    n.chains = 3,
                    n.iter = 10000,
                    n.burnin = 3000,
                    n.thin = 10,
                    DIC = TRUE)

Sys.time()

# Check convergence -------------------------------------------------------

mcmcplot(mod$samples)

gelman.diag(mod$samples, multivariate = F)

rhats <- mod$Rhat
# 
parm <- c("b", "b0", "wA", "wB", "var.process", "deviance")
# 
rhat_graph_fun(list =rhats, parms = parm, rhat = 1.1) +
  labs(title = "KNZ LTER bird SAM Rhat")

# Output summary ----------------------------------------------------------

sum <- summary(mod$samples)

saveRDS(sum, here('examples',
                  'data_output',
                  'bird_fundiv',
                  'model_outputs',
                  'bird_fric_betareg_summary_raw.RDS'))



# Q -----------------------------------------------------------------------



# Loop indexing -----------------------------------------------------------

n.data <- nrow(Q)

n.templag <- Q %>%
  dplyr::select(TAVE:TAVE_l3) %>%
  ncol()

n.npplag <- Q %>%
  dplyr::select(NPP:NPP_l3) %>%
  ncol()

n.transects <- length(unique(Q$TransID))


# Random effects ----------------------------------------------------------

Transect.ID <- Q %>%
  dplyr::select(TransID) %>%
  as_vector()

# Response data -----------------------------------------------------------

#transform like: https://stats.stackexchange.com/questions/31300/dealing-with-0-1-values-in-a-beta-regression
#bray <- as.vector((all_data$observed_all*(n.data-1) + 0.5)/n.data)
fd <- as.vector((Q$Q*(n.data-1) + 0.5)/n.data)

# Covariates --------------------------------------------------------------

Temp <- Q %>%
  dplyr::select(yrID, TransID, TAVE:TAVE_l3) %>%
  pivot_longer(TAVE:TAVE_l3,
               names_to = 'lag',
               values_to = 'temp') %>%
  mutate(temp = scale(temp)) %>%
  pivot_wider(names_from = 'lag',
              values_from = "temp") %>%
  dplyr::select(TAVE:TAVE_l3) %>%
  as.matrix()

sum(is.na(Temp))/(sum(is.na(Temp)) + sum(!is.na(Temp)))

PPT <- Q %>%
  dplyr::select(yrID, TransID, PPT:PPT_l3) %>%
  pivot_longer(PPT:PPT_l3,
               names_to = 'lag',
               values_to = 'ppt') %>%
  mutate(ppt = scale(ppt)) %>%
  pivot_wider(names_from = 'lag',
              values_from = "ppt") %>%
  dplyr::select(PPT:PPT_l3) %>%
  as.matrix()

sum(is.na(PPT))/(sum(is.na(PPT)) + sum(!is.na(PPT)))

# 0% missing data for PPT and TEmp

NPP <- Q %>%
  dplyr::select(yrID, TransID, NPP:NPP_l3) %>%
  pivot_longer(NPP:NPP_l3,
               names_to = 'lag',
               values_to = 'NPP') %>%
  mutate(NPP = scale(NPP)) %>%
  pivot_wider(names_from = 'lag',
              values_from = "NPP") %>%
  dplyr::select(NPP:NPP_l3) %>%
  as.matrix()

sum(is.na(NPP))/(sum(is.na(NPP)) + sum(!is.na(NPP)))

#with subset data, 10% missing

# Make data list ----------------------------------------------------------

data2 <- list(n.data = n.data,
             n.templag = n.templag,
             n.npplag = n.npplag,
             n.transects = n.transects,
             fd = fd,
             Transect.ID = Transect.ID,
             Temp = Temp,
             PPT = PPT,
             NPP = NPP)


Sys.time()
mod2 <- jagsUI::jags(data = data2,
                     inits = NULL,
                     model.file = model,
                     parameters.to.save = params,
                     parallel = TRUE,
                     n.chains = 3,
                     n.iter = 7000,
                     n.burnin = 3000,
                     DIC = TRUE)

Sys.time()

# Check convergence -------------------------------------------------------

mcmcplot(mod2$samples)

gelman.diag(mod2$samples, multivariate = F)

rhats <- mod2$Rhat

parm <- c("b", "b0", "wA", "wB", "var.process", "deviance")

rhat_graph_fun(list =rhats, parms = parm, rhat = 1.1) +
  labs(title = "KNZ LTER bird SAM Rhat")

# Output summary ----------------------------------------------------------

sum <- summary(mod2$samples)

saveRDS(sum, here('examples',
                  'data_output',
                  'bird_fundiv',
                  'model_outputs',
                  'bird_q_betareg_summary_raw.RDS'))







