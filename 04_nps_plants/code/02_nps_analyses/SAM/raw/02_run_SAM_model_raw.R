
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

source(here('00_functions',
            'plot_functions.R'))

# Load data ---------------------------------------------------------------

data_list <- readRDS(here('04_nps_plants',
                          "data_outputs",
                          'SAM',
                          "model_inputs",
                          "nps_diss_SAM_input_data_raw.RDS"))


# Parameters to save ------------------------------------------------------

params <- c('b0.quad',
            'b0.transect',
            #'b0.plot',
            'b',
            'wA',
            'wB',
            'sig.quad',
            'sig.transect',
            #'sig.plot',
            'var.process')



# JAGS model --------------------------------------------------------------

model <- here('04_nps_plants',
              "code", 
              "02_nps_analyses",
              'SAM',
              'raw',
              "jags",
              "plant_SAM_old_raw.R")

Sys.time()
mod <- jagsUI::jags(data = data_list,
                    inits = NULL,
                    model.file = model,
                    parameters.to.save = params,
                    parallel = TRUE,
                    n.chains = 3,
                    n.burnin = 15000,
                    n.thin = 2, 
                    n.iter = 25000, #4000,
                    DIC = TRUE)

Sys.time()

# # Check convergence -------------------------------------------------------
# 
mcmcplot(mod$samples)
# 
gelman.diag(mod$samples, multivariate = F)
# 
rhat_graph_fun(list = mod$Rhat, parms = params, rhat = 1.1)
