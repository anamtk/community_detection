# Running the MSAM for birds
# Ana Miller-ter Kuile
# September 8, 2023

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

data <- readRDS(here("data_outputs",
                     "konza_birds",
                     "model_inputs",
                     "bird_msam_dynmultisite.RDS"))

# Parameters to save ------------------------------------------------------

params <- c(
            #COMMUNITY parameters
            'a1.Effort',
            'a2.Size',
            'lambda.mean',
            'sig.lambda',
            'a0.mean',
            'sig.a0')


#we found ymax to set initials, since otherwise the model will hate us
inits <- function() list(N = data$ymax)#,
                         #omega = data$omega.init)

# JAGS model --------------------------------------------------------------

model <- here("code", 
              "konza_birds",
              "02_kb_analyses",
              'jags',
              "kb_dyn_MSAM_cov.R")

Sys.time()
mod <- jagsUI::jags(data = data,
                         inits = inits,
                         model.file = model,
                         parameters.to.save = params,
                         parallel = TRUE,
                         n.chains = 3,
                         n.iter = 1,
                         DIC = TRUE)

Sys.time()
