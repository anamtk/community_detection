# Running the MSOM for birds
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

data <- readRDS(here('04_nps_plants',
                     'data_outputs',
                     'MSAM',
                     'model_inputs',
                     'nps_msam_dynmultisite.RDS'))

data <- list(n.species = n.species,
             n.quads = n.quads,
             n.yr = n.yr,
             n.rep = n.rep,
             y = y,
             z = z,
             R = R,
             omega.init1 = omega.init1,
             omega.init2 = omega.init2,
             omega.init3 = omega.init3)


# Set Initials ------------------------------------------------------------

inits <- list(list(N = data$z,
                   omega = data$omega.init1),
              list(N = data$z,
                   omega = data$omega.init2),
              list(N = data$z,
                   omega = data$omega.init3))

#ideally you would use the omega.init object to set inits, but it 
#seems to break something with this dataset in particular,
#let me know if this is happening to you and we can troubleshoot
#inits <- function() list(omega = omega.init)

# Parameters to save ------------------------------------------------------

# params <- c(
#   'a1.Effort',
#   'lambda.mean',
#   'sig.lambda',
#   'a0.mean',
#   'sig.a0')

params <- c(
  'psi.mean',
  'sig.lpsi',
  'p.mean',
  'sig.lp',
  'omega'
)


# JAGS model --------------------------------------------------------------

model <- here("04_nps_plants",
              'code',
              '02_nps_analyses',
              'jags',
              "nps_dyn_MSOM_cov.R")

Sys.time()
start<-proc.time()
mod <- jagsUI::jags(data = data,
                    inits = inits,
                    #inits = NULL,
                    model.file = model,
                    parameters.to.save = params,
                    parallel = TRUE,
                    n.chains = 3,
                    n.iter = 1,
                    DIC = TRUE)

Sys.time()
end<-proc.time()
end-start
