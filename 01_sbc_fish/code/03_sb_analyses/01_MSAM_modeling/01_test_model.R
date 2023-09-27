# Running the MSAM for fish
# Ana Miller-ter Kuile
# July 25, 2023

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

data <- readRDS(here('01_sbc_fish',
                     "data_outputs",
                     'MSAM',
                     "model_inputs",
                     "fish_msam_dynmultisite.RDS"))

# Parameters to save ------------------------------------------------------

params <- c(
            #COMMUNITY parameters
            'a1.Vis',
            'a2.Size',
            'lambda.mean',
            'sig.llambda',
            'a0.mean',
            'sig.a0',
            'omega')

#params <- c("N")

#we found ymax to set initials, since otherwise the model will hate us
inits <- list(list(N = data$ymax,
                   omega = data$omega.init),
              list(N = data$ymax,
                   omega = data$omega.init),
              list(N = data$ymax,
                   omega = data$omega.init))

#inits <- function()list(N = data$ymax)

# JAGS model --------------------------------------------------------------

model <- here("01_sbc_fish",
              "code", 
              "03_sb_analyses",
              '01_MSAM_modeling',
              "models",
              "dyn_MSAM_multisite_cov.R")

st.time <- Sys.time()
mod <- jagsUI::jags(data = data,
                    inits = inits,
                    model.file = model,
                    parameters.to.save = params,
                    parallel = TRUE,
                    n.chains = 3,
                    n.iter = 1,
                    DIC = TRUE)

end.time <- Sys.time()

end.time - st.time


# Model without covariance ------------------------------------------------


# Parameters to save ------------------------------------------------------

params <- c(
  #COMMUNITY parameters
  'a1.Vis',
  'a2.Size',
  'lambda.mean',
  'sig.llambda',
  'a0.mean',
  'sig.a0',
  'sig.lambda')


#we found ymax to set initials, since otherwise the model will hate us
# inits <- list(list(N = data$ymax),
#               list(N = data$ymax),
#               list(N = data$ymax))

inits <- function()list(N = data$ymax)

# JAGS model --------------------------------------------------------------

model <- here("01_sbc_fish",
              "code", 
              "03_sb_analyses",
              '01_MSAM_modeling',
              "models",
              "dyn_MSAM_multisite.R")

st.time <- Sys.time()
mod <- jagsUI::jags(data = data,
                    inits = inits,
                    model.file = model,
                    parameters.to.save = params,
                    parallel = TRUE,
                    n.chains = 3,
                    n.iter = 1,
                    DIC = TRUE)

end.time <- Sys.time()

end.time - st.time








