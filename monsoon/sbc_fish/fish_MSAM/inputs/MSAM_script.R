#Monsoon script - fish MSAM
# Ana Miller-ter Kuile
# July 27, 2023

#this script runs the fish MSOM

# Load packages ---------------------------------------------------------------
Sys.time()


# Load packages,
package.list <- c("jagsUI", "coda",
                  'dplyr', 
                  'magrittr', 'tidyr',
                  'mcmcplots') 


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Load Data ---------------------------------------------------------------

#load the formatted data for the JAGS model
data <- readRDS("/scratch/atm234/sbc_fish/inputs/fish_msam_dynmultisite.RDS")

# Compile data ------------------------------------------------------------
data_list <- list(y = data$y,
             vis = data$vis,
             size = data$size,
             n.species = data$n.species,
             n.years = data$n.years,
             n.start = data$n.start,
             n.end = data$n.end,
             n.transects = data$n.transects,
             n.rep = data$n.rep,
             #for initials
             ymax = data$ymax,
             #for covariance prior
             R = data$R)

# Parameters to save ------------------------------------------------------

params <- c(
  #COMMUNITY parameters
  'a1.Vis',
  'a2.Size',
  'lambda.mean',
  'sig.lambda',
  'a0.mean',
  'sig.a0')

# INits -------------------------------------------------------------------

#we found ymax to set initials, since otherwise the model will hate us
inits <- function() list(N = data$ymax)

# JAGS model --------------------------------------------------------------

mod <- jagsUI::jags(data = data_list,
                        inits = inits,
                        #inits = NULL,
                        model.file = '/scratch/atm234/sbc_fish/inputs/dyn_MSAM_multisite_cov.R',
                        parameters.to.save = params,
                        parallel = TRUE,
                        n.chains = 3,
                        n.iter = 4000,
                        n.burnin = 100,
                        DIC = TRUE)

#save as an R data object
saveRDS(mod, 
        file ="/scratch/atm234/sbc_fish/outputs/fish_MSAM_model.RDS")

Sys.time()



# Check convergence -------------------------------------------------------

mcmcplot(mod$samples,
         dir = "/scratch/atm234/sbc_fish/outputs/mcmcplots/MSAM")

# Get RHat per parameter ------------------------------------------------

Rhat <- mod$Rhat

saveRDS(Rhat, "/scratch/atm234/sbc_fish/outputs/fish_MSAM_model_Rhat.RDS")


# Get Raftery diag --------------------------------------------------------


raf <- raftery.diag(mod$samples)

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



