
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
            'wB',
            'wC')


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


# Temp --------------------------------------------------------------------


# JAGS model --------------------------------------------------------------

model <- here("code", 
              "03_analyses",
              '03_metric_SAM_modeling',
              "jags",
              "gain_SAM_temp.R")

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


# chla --------------------------------------------------------------------


# JAGS model --------------------------------------------------------------

model <- here("code", 
              "03_analyses",
              '03_metric_SAM_modeling',
              "jags",
              "gain_SAM_chla.R")

Sys.time()
modg2 <- jagsUI::jags(data = data_list,
                     inits = NULL,
                     model.file = model,
                     parameters.to.save = params,
                     parallel = TRUE,
                     n.chains = 3,
                     n.iter = 4000,
                     DIC = TRUE)

Sys.time()

# Check convergence -------------------------------------------------------

mcmcplot(modg2$samples)

# Model selection ---------------------------------------------------------

parms <- c("lpd", "pd", "Dsum")

t_selection <- update(modg,
                      parameters.to.save = parms,
                      n.iter = 500)

c_selection <- update(modg2,
                      parameters.to.save = parms,
                      n.iter= 500)


# WAIC --------------------------------------------------------------------

t_sum <- summary(t_selection$samples)
c_sum <- summary(c_selection$samples)

#get pd rows
lppd_t <- sum(log(t_sum$statistics[663:1324, 1]))
lppd_c <- sum(log(c_sum$statistics[663:1324,1]))

#get lpd columsn
pwaic_t <- sum(t_sum$statistics[1:662, 2]^2)
pwaic_c <- sum(c_sum$statistics[1:662, 2]^2)

(waic_t <- lppd_t + pwaic_t)
(waic_c <- lppd_c + pwaic_c)

# DIC ---------------------------------------------------------------------

modg$DIC
modg2$DIC

# Dinf --------------------------------------------------------------------

t_selection$mean$Dsum
c_selection$mean$Dsum

# Look at some plots ------------------------------------------------------

# sumg <- summary(modg$samples)
# 
# medg <- as.data.frame(sumg$quantiles) %>%
#   rownames_to_column(var = "parameter") %>%
#   dplyr::select(parameter, `2.5%`, `50%`, `97.5%`) %>%
#   filter(parameter != "deviance")
# 
# write.csv(medg, here("data_outputs",
#                     "SAM_outputs",
#                     "gain_SAM_summary.csv"))




