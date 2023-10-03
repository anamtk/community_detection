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
            'S',
            'G')

#params <- c("N")

#we found ymax to set initials, since otherwise the model will hate us
# inits <- list(list(N = data$ymax,
#                    omega = data$omega.init1),
#               list(N = data$ymax,
#                    omega = data$omega.init2),
#               list(N = data$ymax,
#                    omega = data$omega.init3))

inits <- list(list(N = data$ymax),
              list(N = data$ymax),
              list(N = data$ymax))


#inits <- function()list(N = data$ymax)

# JAGS model --------------------------------------------------------------

model <- here("01_sbc_fish",
              "code", 
              "03_sb_analyses",
              '01_MSAM_modeling',
              "models",
              "dyn_MSAM_multisite.R")

(st.time <- Sys.time())
mod <- jagsUI::jags(data = data,
                    #inits = NULL,
                    inits = inits,
                    model.file = model,
                    parameters.to.save = params,
                    parallel = TRUE,
                    n.chains = 3,
                    n.iter = 4000,
                    DIC = TRUE)

end.time <- Sys.time()

end.time - st.time

mcmcplot(mod$samples)

# Get inits ---------------------------------------------------------------

#get the MCMC chains
samples <- mod$samples

#function to make each chain a dataframe
df_fun <- function(chain){
  df <- as.data.frame(chain) %>%
    rownames_to_column(var = "iteration")
  return(df)
}

#use that function on all list elements
samp_dfs <- lapply(samples, df_fun)

#make into one dataframe
samp_df <- bind_rows(samp_dfs, .id = "chain")

#get values for all parameters from the last iteration of the
#chain with the lowest deviance
samp_df2 <- samp_df %>%
  group_by(chain) %>%
  #get mean deviance by chain
  mutate(mean_dev = mean(deviance, na.rm = T)) %>%
  ungroup() %>%
  #get only the chain with the minimum average deviance
  filter(mean_dev == min(mean_dev)) %>%
  #pull out the final iteration from that chain
  filter(iteration == max(iteration)) %>%
  dplyr::select(-chain, -iteration,
                -deviance, -mean_dev) 

#for root nodes (except species-level taus)
a1.Vis <- as.vector(samp_df2$a1.Vis)
a2.Size <- as.vector(samp_df2$a2.Size)
lambda.mean <- as.vector(samp_df2$lambda.mean)
sig.llambda <- as.vector(samp_df2$sig.llambda)
a0.mean <- as.vector(samp_df2$a0.mean)
sig.a0 <- as.vector(samp_df2$sig.a0)


# Run model again ---------------------------------------------------------

inits2 <- list(list(N = data$ymax,
                   a1.Vis = a1.Vis,
                   a2.Size = a2.Size,
                   lambda.mean = lambda.mean,
                   sig.llambda = sig.llambda,
                   a0.mean = a0.mean,
                   sig.a0 = sig.a0),
              list(N = data$ymax,
                   a1.Vis = a1.Vis + 0.5,
                   a2.Size = a2.Size + 0.5,
                   lambda.mean = lambda.mean + 0.1,
                   sig.llambda = sig.llambda + 0.001,
                   a0.mean = a0.mean + 0.001,
                   sig.a0 = sig.a0 + 5),
              list(N = data$ymax,
                   a1.Vis = a1.Vis - 0.5,
                   a2.Size = a2.Size - 0.5,
                   lambda.mean = lambda.mean - 0.1,
                   sig.llambda = sig.llambda - 0.001,
                   a0.mean = a0.mean - 0.001,
                   sig.a0 = sig.a0 - 5))


# Run again ---------------------------------------------------------------


(st.time <- Sys.time())
mod2 <- jagsUI::jags(data = data,
                    inits = inits2,
                    model.file = model,
                    parameters.to.save = params,
                    parallel = TRUE,
                    n.chains = 3,
                    n.iter = 4000,
                    DIC = TRUE)

end.time <- Sys.time()

end.time - st.time

mcmcplot(mod2$samples)

