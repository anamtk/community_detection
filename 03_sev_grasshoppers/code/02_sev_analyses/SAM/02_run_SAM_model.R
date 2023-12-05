
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


# Load data ---------------------------------------------------------------

data_list <- readRDS(here('03_sev_grasshoppers',
                          "data_outputs",
                          'SAM',
                          "model_inputs",
                          "sev_bray_SAM_input_data.RDS"))


# Parameters to save ------------------------------------------------------

params <- c('b0.web',
            'b0.transect',
            'b',
            'wA',
            'wB',
            'wC',
            'sig.web',
            'sig.transect',
            'var.process')



# JAGS model --------------------------------------------------------------

model <- here('03_sev_grasshoppers',
              "code", 
              "02_sev_analyses",
              'SAM',
              "jags",
              "sev_SAM_2.R")

Sys.time()
mod <- jagsUI::jags(data = data_list,
                    inits = NULL,
                    model.file = model,
                    parameters.to.save = params,
                    parallel = TRUE,
                    n.chains = 3,
                    n.iter = 1,
                    DIC = TRUE)

Sys.time()

# # Check convergence -------------------------------------------------------
# 
mcmcplot(mod$samples)
# 
gelman.diag(mod$samples, multivariate = F)
# 
#
