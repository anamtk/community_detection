
# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 
package.list <- c("here", "tidyverse", 
                  "jagsUI",
                  'rjags',
                  'mcmcplots',
                  "coda",
                  'patchwork') #mcmc output


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Gains -------------------------------------------------------------------

# Load data ---------------------------------------------------------------

data_list <- readRDS(here('04_nps_plants',
                          "data_outputs",
                          'SAM',
                          "model_inputs",
                          "nps_gain_SAM_input_data.RDS"))

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

model <- here('04_nps_plants',
              "code", 
              "02_nps_analyses",
              'SAM',
              'modeled',
              "jags",
              "plant_SAM_old.R")

Sys.time()
mod <- jagsUI::jags(data = data_list,
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

saveRDS(sum, here('04_nps_plants',
                  'data_outputs',
                  'SAM',
                  'model_outputs',
                  'nps_gain_SAM_summary.RDS'))



# Losses ------------------------------------------------------------------

# Load data ---------------------------------------------------------------

data_list2 <- readRDS(here('04_nps_plants',
                          "data_outputs",
                          'SAM',
                          "model_inputs",
                          "nps_loss_SAM_input_data.RDS"))


# JAGS model --------------------------------------------------------------

Sys.time()
mod2 <- jagsUI::jags(data = data_list2,
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

saveRDS(sum2, here('04_nps_plants',
                  'data_outputs',
                  'SAM',
                  'model_outputs',
                  'nps_loss_SAM_summary.RDS'))
