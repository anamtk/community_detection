
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
                          "bray_SAM_input_data_long.RDS"))


# Parameters to save ------------------------------------------------------

params <- c('b0.transect',
            'b',
            'wA',
            'wB',
            'diff',
            'sig.transect',
            'sig.site',
            'var.process')



# JAGS model --------------------------------------------------------------

model <- here('01_sbc_fish',
              "code", 
              "03_sb_analyses",
              '02_SAM_modeling',
              'modeled',
              "jags",
              "fish_SAM.R")

Sys.time()
mod <- jagsUI::jags(data = data_list,
                    inits = NULL,
                    model.file = model,
                    parameters.to.save = params,
                    parallel = TRUE,
                    n.chains = 3,
                    n.burnin =1500,
                    n.iter = 5000,
                    DIC = TRUE)

Sys.time()

# Check convergence -------------------------------------------------------
# 
mcmcplot(mod$samples)
# 
gelman.diag(mod$samples, multivariate = F)

rhat <- mod$Rhat

parm <- c('b0.transect', 'b', 'wA', 'wB', 'sig.transect', 'sig.site',
          'var.process', 'deviance')
rhat_graph_fun(rhat, parms = parm, rhat = 1.1) +
  labs(title = "SBC LTER fish SAM Rhat")

ggsave(plot = last_plot(),
       filename = here("pictures",
                       "supplementary",
                       'SAM',
                       "SBC_SAM_Rhat_graph.jpg"),
       height = 4,
       width = 6,
       units = "in")
# 
# # Look at some plots ------------------------------------------------------
# 
sum <- summary(mod$samples)
# 
saveRDS(sum, here("01_sbc_fish",
                  "data_outputs",
                  "SAM",
                  "model_outputs",
                  "fish_SAM_summary.RDS"))
# 

# GOF ---------------------------------------------------------------------

parms2 <- c("bray.rep", "resid")

mod2 <- update(mod,
               parameters.to.save = parms2,
               n.iter = 500)

sum2 <- summary(mod2$samples)

saveRDS(sum2, here("01_sbc_fish",
                   "data_outputs",
                   "SAM",
                   "model_outputs",
                   "fish_SAM_GOF_summary.RDS"))
