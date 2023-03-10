# Running the nest initiation
# Ana Miller-ter Kuile
# November 4, 2021

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

# Fix parallels -----------------------------------------------------------

#hopefully the parallels issue gets fixed, but for now this if statement works
# to set system preferences for jags to run with parallel
if (Sys.getenv("RSTUDIO") == "1" && !nzchar(Sys.getenv("RSTUDIO_TERM")) && 
    Sys.info()["sysname"] == "Darwin" && getRversion() >= "4.0.0") {
  parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
}


# Load Data ---------------------------------------------------------------

data <- readRDS(here("data_outputs",
                     "model_inputs",
                     "fish_data_singlesite.RDS"))

# Parameters to save ------------------------------------------------------

params <- c(#species-level parameters
            "eta",
            "a1.Vis",
            'a2.Size',
            #community-level parameters
            'tau.eta',
            'sd.vis',
            'sd.size',
            'sd.lpsi',
            'sd.lp',
            "rho")


# JAGS model --------------------------------------------------------------

model <- here("code", 
              "03_analyses",
              'occupancy_modeling',
              "models",
              "test_model.R")

#this model ran in 15 minutes on my desktop,
#but will likely take longer when we add more sites
Sys.time()

jags <- jagsUI::jags(data = data,
                         inits = NULL,
                         model.file = model,
                         parameters.to.save = params,
                         parallel = TRUE,
                         n.chains = 3,
                         n.iter = 60000,
                         n.burnin = 1000,
                         DIC = TRUE)

Sys.time()

# Check convergence -------------------------------------------------------

mcmcplot(jags$samples)


# Raftery -----------------------------------------------------------------

raf <- raftery.diag(jags$samples)

names <- rownames(raf[[1]]$resmatrix)
ch1 <- raf[[1]]$resmatrix[,2]
ch2 <- raf[[2]]$resmatrix[,2]
ch3 <- raf[[3]]$resmatrix[,2]

raf_all <- as.data.frame(cbind(names, 
                               ch1, ch2, ch3)) %>%
  mutate(ch1 = as.numeric(ch1),
         ch2 = as.numeric(ch2),
         ch3 = as.numeric(ch3)) %>%
  filter(!str_detect(names, "z")) %>%
  pivot_longer(ch1:ch3,
               names_to = "chain",
               values_to = 'iterations') 

ggplot(raf_all, aes(x = iterations/3)) +
  geom_histogram() 

raf_all %>%
  summarise(iterations_90 = quantile(iterations, 
                                     probs = 0.9, 
                                     na.rm = T)/3,
            iterations_95 = quantile(iterations,
                                     probs = 0.95,
                                     na.rm = T)/3,
            max = max(iterations, 
                      na.rm = T)/3)
# A tibble: 1 × 3
# iterations_90 iterations_95   max
# <dbl>         <dbl> <dbl>
#   1        20717.        29211. 86112

bu1 <- raf[[1]]$resmatrix[,1]
bu2 <- raf[[2]]$resmatrix[,1]
bu3 <- raf[[3]]$resmatrix[,1]

burn <- as.data.frame(cbind(names, bu1, bu2, bu3)) %>%
  mutate(bu1 = as.numeric(bu1),
         bu2 = as.numeric(bu2),
         bu3 = as.numeric(bu3)) %>%
  filter(!str_detect(names, "z")) %>%
  pivot_longer(bu1:bu3,
               names_to = "chain",
               values_to = 'iterations') 

burn %>%
  summarise(max(iterations, na.rm = T))
#792


# Long run outputs --------------------------------------------------------

#did they converge? only tau.eta looks funky here
Rhat <- jags$Rhat

#get a summary of all the output parameters that we tracked
summary <- jags$summary

#update converged model to get psi values
psi <- update(jags,
              parameters.to.save = c("psi"),
              n.iter = 5000)

#summary of psi run
psi_sum <- psi$summary

community_params <- update(jags,
                           parameters.to.save = c("mu.vis",
                                                  'sd.vis',
                                                  'mu.size',
                                                  'sd.size'),
                           n.iter = 5000)

community_sum <- community_params$summary


output_summaries <- list(Rhat = Rhat,
                         summary = summary,
                         psi_sum = psi_sum,
                         community_sum = community_sum)



saveRDS(output_summaries, here("data_outputs",
                               "monsoon_outputs",
                               "fish_MSOM_1_26_stats.RDS"))


# Export the corrected 1-0 matrices ---------------------------------------

#update converged model to get z-matrix values
z_matrix <- update(jags,
                   parameters.to.save = c("z"),
                   n.iter = 500)

saveRDS(z_matrix$samples, here("data_outputs",
                       "community_matrices",
                       "fish_MSOM_community_matrices.RDS"))

#get summary stats for that
# z_50 <- z_matrix$q50$z
# z_2.5 <- z_matrix$q2.5$z
# z_97.5 <- z_matrix$q97.5$z
# 
# z_matrices <- list(z_50 = z_50,
#                    z_2.5 = z_2.5,
#                    z_97.5 = z_97.5)
# 
# saveRDS(z_matrices, here("data_outputs",
#                          "monsoon_outputs",
#                          "fish_MSOM_1_26_matrices.RDS"))
