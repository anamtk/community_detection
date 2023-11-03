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

data <- readRDS(here("02_konza_birds",
                     "data_outputs",
                     "MSAM",
                     "model_inputs",
                     "bird_msam_dynmultisite.RDS"))

# Parameters to save ------------------------------------------------------

params <- c(
            #COMMUNITY parameters
            'a1.Effort',
            'a2.Size',
            'mu.llambda',
            'sig.llambda',
            'mu.a0',
            'sig.a0')


#we found ymax to set initials, since otherwise the model will hate us
inits <- function() list(N = data$ymax)#,
                         #omega = data$omega.init)

# JAGS model --------------------------------------------------------------

model <- here("02_konza_birds",
              "code", 
              "02_kb_analyses",
              "MSAM",
              'jags',
              "kb_MSAM_simple_noeffort.R")

start.time <- Sys.time()
mod <- jagsUI::jags(data = data,
                         inits = inits,
                         model.file = model,
                         parameters.to.save = params,
                         parallel = TRUE,
                         n.chains = 3,
                         n.iter = 4000,
                         DIC = TRUE)

end.time <- Sys.time()

end.time - start.time

mcmcplot(mod$samples)
gelman.diag(mod$samples)

mod2 <- update(mod,
               parallel = TRUE,
               n.iter = 4000)

mcmcplot(mod2$samples)
gelman.diag(mod2$samples)

mod3 <- update(mod2,
               parallel = TRUE,
               n.iter = 4000)

mcmcplot(mod3$samples)
gelman.diag(mod3$samples)

mod4 <- update(mod3, 
               parallel = TRUE,
               n.iter = 4000)

mcmcplot(mod4$samples)
gelman.diag(mod4$samples)
