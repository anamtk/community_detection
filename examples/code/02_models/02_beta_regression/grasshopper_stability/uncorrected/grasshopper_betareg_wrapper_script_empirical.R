
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

source(here('examples',
            'code',
            '00_functions',
            'plot_functions.R'))
# Load data ---------------------------------------------------------------

data_list <- readRDS(here('examples',
                          'data_output',
                          'grasshopper_stability',
                          'model_inputs',
                          'grasshopper_betareg_input_data_list_raw.RDS'))


# Parameters to save ------------------------------------------------------

params <- c('b0.web',
            'b0.transect',
            'b',
            'b0',
            'wA',
            'wB',
            'wC',
            'deltaA',
            'deltaB',
            'deltaC',
            'sig.web',
            'sig.transect',
            'var.process',
            'cumm.tempwt',
            'cumm.pptwt',
            'cumm.nppwt')



# JAGS model --------------------------------------------------------------

model <- here('examples',
              'code',
              '02_models',
              '02_beta_regression',
              'grasshopper_stability',
              'uncorrected',
              'grasshopper_betareg_model_empirical.R')

Sys.time()
mod <- jagsUI::jags(data = data_list,
                    inits = NULL,
                    model.file = model,
                    parameters.to.save = params,
                    parallel = TRUE,
                    n.chains = 3,
                    n.burnin = 2500,
                    n.iter = 5000,
                    DIC = TRUE)

Sys.time()

# # Check convergence -------------------------------------------------------
# 
mcmcplot(mod$samples)
# 
gelman.diag(mod$samples, multivariate = F)
# 
#
rhat <- mod$Rhat

parms <- c('b0.web', 'b0.transect',
           'b', 'b0', 'wA', 'wB', "wC",
           'sig.web', 'sig.transect', 'var.process',
           'deviance')

rhat_graph_fun(rhat, parms = parms, rhat = 1.1) +
  labs(title = "SEV LTER grasshopper SAM Rhat: \n raw data")



# Summary -----------------------------------------------------------------

sum <- summary(mod$samples)

saveRDS(sum(here('examples',
                 'data_output',
                 'grasshopper_stability',
                 'model_outputs',
                 'grasshopper_betareg_model_summary_raw.RDS')))
