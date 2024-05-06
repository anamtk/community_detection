
#Monsoon script - nps_plants
# November 17, 2023

#this script runs the grasshopper SAM

# Load packages ---------------------------------------------------------------
(start.time <- Sys.time())


# Load packages
package.list <- c("jagsUI", "coda",
                  'dplyr', 'stringr',
                  'magrittr', 'tidyr',
                  'mcmcplots','ggplot2') 


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

# Load Data ---------------------------------------------------------------

#load the formatted data for the JAGS model
data <- readRDS("/scratch/sml665/nps_plants/SAM/inputs/nps_diss_SAM_input_data.RDS")


# Compile data ------------------------------------------------------------

data_list <- list(n.data = data$n.data,
                  n.lag = data$n.lag,
                  n.quads = data$n.quads,
                  n.transects = data$n.transects,
                  n.plots = data$n.plots,
                  diss = data$diss,
                  var.estimate = data$var.estimate,
                  Plot.ID = data$Plot.ID,
                  Transect.ID = data$Transect.ID,
                  Quad.ID = data$Quad.ID,
                  PPT = data$PPT, 
                  VPD = data$VPD)


# Parameters to save ------------------------------------------------------

params <- c(#'b0.quad',
  'b0.transect',
  #'b0.plot',
  'b',
  'wA',
  'wB',
  #'sig.quad',
  'sig.transect',
  #'sig.plot',
  'var.process')

# INits -------------------------------------------------------------------

#inits <- readRDS("/scratch/atm234/whwo_ipm/parameter_models/adult_occupancy/inputs/adult_occupancy_inits.RDS")

# JAGS model --------------------------------------------------------------

mod <- jagsUI::jags(data = data_list,
                    #inits = inits,
                    inits = NULL,
                    model.file = '/scratch/sml665/nps_plants/SAM/inputs/plant_SAM_old.R',
                    parameters.to.save = params,
                    parallel = TRUE,
                    n.chains = 3,
                    n.burnin = 1000,
                    n.iter =  41000,
                    n.thin = 7,
                    DIC = TRUE)

#save as an R data object
saveRDS(mod, 
        file ="/scratch/sml665/nps_plants/SAM/outputs/nps_SAM_model.RDS")

(end.time <- Sys.time())

(tot.time <- end.time - start.time)

# Check convergence -------------------------------------------------------

mcmcplot(mod$samples,
         dir = "/scratch/sml665/nps_plants/SAM/outputs/mcmcplots/SAM")

# Get RHat per parameter ------------------------------------------------

Rhat <- mod$Rhat

saveRDS(Rhat, "/scratch/sml665/nps_plants/SAM/outputs/nps_SAM_model_Rhat.RDS")


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
  pivot_longer(ch1:ch3,
               names_to = "chain",
               values_to = 'iterations') 

raf_all %>%
  summarise(iterations_90 = quantile(iterations, 
                                     probs = 0.9, 
                                     na.rm = T)/3,
            iterations_95 = quantile(iterations,
                                     probs = 0.95,
                                     na.rm = T)/3,
            max = max(iterations, 
                      na.rm = T)/3)

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
