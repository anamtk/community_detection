
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

data_list <- readRDS(here('01_sbc_fish',
                          "data_outputs",
                          'SAM',
                          "model_inputs",
                          "bray_SAM_input_data.RDS"))


# Parameters to save ------------------------------------------------------

params <- c('b0.transect',
            'b',
            'wA',
            'wB',
            'wC',
            'sig.transect',
            'sig.site',
            'var.process')



# JAGS model --------------------------------------------------------------

model <- here('01_sbc_fish',
              "code", 
              "03_sb_analyses",
              '02_SAM_modeling',
              "jags",
              "fish_SAM.R")

Sys.time()
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

# alpha and beta explorations ---------------------------------------------

alphaX <- as.data.frame(mod$sims.list$alphaX) %>%
  pivot_longer(cols = everything(),
               names_to = "i",
               values_to = "alphaX")

aX <- alphaX %>%
  filter(alphaX < 1) %>%
  tally() %>%
  as_vector()

aX/7944000
#55% of the time, alphaX is <1

betaX <- as.data.frame(mod$sims.list$betaX) %>%
  pivot_longer(cols = everything(),
               names_to = "i",
               values_to = "betaX")

bX <- betaX %>%
  filter(betaX < 1) %>%
  tally() %>%
  as_vector()

bX/7944000
#57% of the time, betaX <1

ggplot(betaX, aes(x = betaX)) +
  geom_histogram() +
  geom_vline(xintercept = 1)

ggplot(alphaX, aes(x = alphaX)) +
  geom_histogram() +
  geom_vline(xintercept = 1)

diff <- as.data.frame(mod$sims.list$diff) %>%
  pivot_longer(cols = everything(),
               names_to = "i",
               values_to = "diff")

diffX <- diff %>%
  filter(diff < 0) %>%
  tally() %>%
  as_vector()

diffX/7944000

ggplot(diff, aes(x = diff)) +
  geom_histogram() +
  geom_vline(xintercept = 0)

# Look at some plots ------------------------------------------------------

sum <- summary(mod$samples)

med <- as.data.frame(sum$quantiles) %>%
  rownames_to_column(var = "parameter") %>%
  dplyr::select(parameter, `2.5%`, `50%`, `97.5%`) %>%
  filter(parameter != "deviance")

write.csv(med, here("data_outputs",
                    "SAM_outputs",
                    "loss_SAM_summary.csv"))

