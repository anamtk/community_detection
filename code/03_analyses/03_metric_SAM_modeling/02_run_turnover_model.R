
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

data_list <- readRDS(here("data_outputs",
                          "model_inputs",
                          "turnover_SAM",
                          "turnover_SAM_input_data.RDS"))


# Parameters to save ------------------------------------------------------

params <- c("b0",
            'b',
            'wA',
            'wB')


# Loss --------------------------------------------------------------------


# JAGS model --------------------------------------------------------------

model <- here("code", 
              "03_analyses",
              '03_metric_SAM_modeling',
              "jags",
              "loss_SAM.R")

Sys.time() #~2 minutes??
mod <- jagsUI::jags(data = data_list,
                    inits = NULL,
                    model.file = model,
                    parameters.to.save = params,
                    parallel = TRUE,
                    n.chains = 3,
                    n.iter = 4000,
                    DIC = TRUE)

Sys.time()

# Check convergence -------------------------------------------------------

mcmcplot(mod$samples)

gelman.diag(mod$samples, multivariate = F)


# Look at some plots ------------------------------------------------------

sum <- summary(mod$samples)

med <- as.data.frame(sum$quantiles) %>%
  rownames_to_column(var = "parameter") %>%
  dplyr::select(parameter, `2.5%`, `50%`, `97.5%`) %>%
  filter(parameter != "deviance")

write.csv(med, here("data_outputs",
                    "SAM_outputs",
                    "loss_SAM_summary.csv"))

# Gain model --------------------------------------------------------------

# JAGS model --------------------------------------------------------------

model <- here("code", 
              "03_analyses",
              '03_metric_SAM_modeling',
              "jags",
              "gain_SAM.R")

Sys.time()
modg <- jagsUI::jags(data = data_list,
                    inits = NULL,
                    model.file = model,
                    parameters.to.save = params,
                    parallel = TRUE,
                    n.chains = 3,
                    n.iter = 4000,
                    DIC = TRUE)

Sys.time()

# Check convergence -------------------------------------------------------

mcmcplot(modg$samples)

gelman.diag(modg$samples, multivariate = F)


# Look at some plots ------------------------------------------------------

sumg <- summary(modg$samples)

medg <- as.data.frame(sumg$quantiles) %>%
  rownames_to_column(var = "parameter") %>%
  dplyr::select(parameter, `2.5%`, `50%`, `97.5%`) %>%
  filter(parameter != "deviance")

write.csv(medg, here("data_outputs",
                    "SAM_outputs",
                    "gain_SAM_summary.csv"))




