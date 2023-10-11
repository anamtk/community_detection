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

n.yr <- as.vector(unname(data$n.yr))

data <- list(n.species = data$n.species,
             n.quads = data$n.quads,
             n.yr = n.yr,
             n.rep = data$n.rep,
             y = data$y,
             z = data$z)


# Set Initials ------------------------------------------------------------

inits <- list(list(N = data$z),
              list(N = data$z),
              list(N = data$z))

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
  'phi.mean',
  'sig.lphi',
  'gamma.mean',
  'sig.gamma',
  'p.mean',
  'sig.lp'
)


# JAGS model --------------------------------------------------------------

model <- here("04_nps_plants",
              'code',
              '02_nps_analyses',
              'jags',
              "nps_dyn_MSOM_simple.R")

Sys.time()
start<-proc.time()
mod <- jagsUI::jags(data = data,
                    inits = inits,
                    model.file = model,
                    parameters.to.save = params,
                    parallel = TRUE,
                    n.chains = 3,
                    n.iter = 1,
                    DIC = TRUE)

Sys.time()
end<-proc.time()
end-start
