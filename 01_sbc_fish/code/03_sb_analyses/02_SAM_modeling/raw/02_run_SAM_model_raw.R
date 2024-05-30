
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

source(here('00_functions',
            'plot_functions.R'))

# Load data ---------------------------------------------------------------

data_list <- readRDS(here('01_sbc_fish',
                          "data_outputs",
                          'SAM',
                          "model_inputs",
                          "bray_SAM_input_data_raw.RDS"))


# Parameters to save ------------------------------------------------------

params <- c('b0.transect',
            'b0.site',
            'b',
            'wA',
            'wB',
            'sig.transect',
            'sig.site',
            'var.process')



# JAGS model --------------------------------------------------------------

model <- here('01_sbc_fish',
              "code", 
              "03_sb_analyses",
              '02_SAM_modeling',
              'raw',
              "jags",
              "raw_fish_SAM.R")

Sys.time()
mod <- jagsUI::jags(data = data_list,
                    inits = NULL,
                    model.file = model,
                    parameters.to.save = params,
                    parallel = TRUE,
                    n.chains = 3,
                    n.iter = 5000,
                    n.burnin = 1500,
                    DIC = TRUE)

Sys.time()

# Check convergence -------------------------------------------------------
# 
mcmcplot(mod$samples)
# 
gelman.diag(mod$samples, multivariate = F)

rhat_graph_fun(list = mod$Rhat, parms = params, rhat = 1.1) +
  labs(title = "SBC LTER fish SAM Rhat: \n observed data")

ggsave(plot = last_plot(),
       filename = here("pictures",
                       "supplementary",
                       'SAM',
                       "SBC_SAM_raw_Rhat_graph.jpg"),
       height = 4,
       width = 6,
       units = "in")
# 

# Save summary ------------------------------------------------------------

sum <- summary(mod$samples)

saveRDS(sum, here("01_sbc_fish",
                  "data_outputs",
                  "SAM",
                  "model_outputs",
                  "fish_SAM_summary_raw.RDS"))


# GOF ---------------------------------------------------------------------

parms2 <- c('bray.rep', 'resid')

mod2 <- update(mod,
               parameters.to.save = parms2,
               n.iter = 500)

sum2 <- summary(mod2$samples)

saveRDS(sum2, here("01_sbc_fish",
                  "data_outputs",
                  "SAM",
                  "model_outputs",
                  "fish_SAM_GOF_summary_raw.RDS"))
