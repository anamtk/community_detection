# Running the MSAM for grasshoppers
# Ana Miller-ter Kuile
# September 11, 2023

# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 
package.list <- c("here", "tidyverse", 
                  "jagsUI",
                  'rjags',
                  'mcmcplots',
                  "coda") #mcmc output


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load Data ---------------------------------------------------------------

data <- readRDS(here('03_sev_grasshoppers',
                     "data_outputs",
                     "MSAM",
                     "model_inputs",
                     "sev_msam_dynmultisite.RDS"))

# Parameters to save ------------------------------------------------------

params <- c(
            #COMMUNITY parameters
            'mu.llambda',
            'sig.llambda',
            'p.mean',
            'sig.lp')


#we found ymax to set initials, since otherwise the model will hate us
inits <- function() list(N = data$ymax)

# JAGS model --------------------------------------------------------------

model <- here('03_sev_grasshoppers',
              "code",
              '02_sev_analyses',
              "MSAM",
              'jags',
              'sev_MSAM_simple_nocov.R')

start.time <- Sys.time()
mod <- jagsUI::jags(data = data,
                         inits = inits,
                         model.file = model,
                         parameters.to.save = params,
                         parallel = TRUE,
                         n.chains = 3,
                         n.iter = 1,
                         DIC = TRUE)

end.time <- Sys.time()

end.time - start.time

mcmcplot(mod$samples)
gelman.diag(mod$samples, multivariate = F)
