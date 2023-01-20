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
                     "data_singlesite_MSOM.RDS"))

# Parameters to save ------------------------------------------------------

params <- c("a1.Vis",
            "sd.vis",
            "eta",
            "tau.eta",
            "rho",
            "lpsi",
            "gamma",
            "mu.vis",
            "sd.vis",
            "mu.gamma",
            "sigmagamma"
            )


# JAGS model --------------------------------------------------------------

model <- here("code", 
              "03_analyses",
              "models",
              "test_dynamic_model.R")

Sys.time()

jags <- jagsUI::jags(data = data,
                         inits = NULL,
                         model.file = model,
                         parameters.to.save = params,
                         parallel = TRUE,
                         n.chains = 3,
                         n.iter = 4000,
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


# Plot visibility by species effect ---------------------------------------

df <- as.data.frame(cbind(
vis.50 = jags$q50$a1.Vis,
vis.lower = jags$q2.5$a1.Vis,
vis.upper = jags$q97.5$a1.Vis,
species = 1:length(jags$q50$a1.Vis)))

df %>%
  filter(species %in% c(1:122)) %>%
ggplot(aes(x = species, y = vis.50)) +
  geom_point() +
  geom_errorbar(aes(ymin = vis.lower, ymax = vis.upper)) +
  coord_flip()



