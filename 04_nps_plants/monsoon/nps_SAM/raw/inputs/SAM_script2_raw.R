#Monsoon script - grasshopper
# Ana Miller-ter Kuile
# May 12, 2023

#this script runs the plant SAM

# Load packages ---------------------------------------------------------------
(start.time <- Sys.time())


# Load packages
package.list <- c("jagsUI", "coda",
                  'dplyr', 'stringr',
                  'magrittr', 'tidyr',
                  'mcmcplots','ggplot2',
                  'tibble', 'purrr') 


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load model --------------------------------------------------------------

mod <- readRDS(file ='/scratch/atm234/nps_plants/SAM/raw/outputs/nps_SAM_model_raw.RDS')

# Get initials from previous model ----------------------------------------

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


var.process <- as.vector(samp_df2$var.process)
sig.quad <- as.vector(samp_df2$sig.quad)
sig.transect <- as.vector(samp_df2$sig.transect)
b0 <- as.vector(samp_df2$b0)
b <- samp_df2 %>%
  dplyr::select(contains('b')) %>%
  dplyr::select(-contains('b0')) %>%
  dplyr::select(-contains('wB')) %>%
  dplyr::select(-contains('delta')) %>%
  dplyr::select(-contains('web')) %>%
  as_vector()

deltaA <- samp_df2 %>%
  dplyr::select(contains('deltaA')) %>%
  as_vector()

deltaB <- samp_df2 %>%
  dplyr::select(contains('deltaB')) %>%
  as_vector()

deltaC <- samp_df2 %>%
  dplyr::select(contains('deltaC')) %>%
  as_vector()



# Load Data ---------------------------------------------------------------

#load the formatted data for the JAGS model
data <- readRDS("/scratch/atm234/nps_plants/SAM/raw/inputs/nps_diss_SAM_input_data_raw.RDS")

# Compile data ------------------------------------------------------------

data_list <- list(n.data = data$n.data,
                  n.lag = data$n.lag,
                  n.quads = data$n.quads,
                  n.transects = data$n.transects,
                  n.plots = data$n.plots,
                  diss = data$diss,
                  Plot.ID = data$Plot.ID,
                  Transect.ID = data$Transect.ID,
                  Quad.ID = data$Quad.ID,
                  PPT = data$PPT, 
                  VPD = data$VPD)

# Parameters to save ------------------------------------------------------

params <- c('b0.quad',
            'b0.transect',
            #'b0.plot',
            'b',
            'wA',
            'wB',
            'sig.quad',
            'sig.transect',
            #'sig.plot',
            'var.process')

# INits -------------------------------------------------------------------

#inits <- readRDS("/scratch/atm234/whwo_ipm/parameter_models/adult_occupancy/inputs/adult_occupancy_inits.RDS")

inits <- list(list(b0 = b0,
                   b = b,
                   sig.quad = sig.quad,
                   sig.transect = sig.transect,
                   deltaA = deltaA,
                   deltaB = deltaB,
                   deltaC = deltaC),
              list(b0 = b0 + 0.05,
                   b = b + 0.5,
                   sig.quad = sig.quad + 0.01,
                   sig.transect = sig.transect + 0.01,
                   deltaA = deltaA + 0.01,
                   deltaB = deltaB + 0.01,
                   deltaC = deltaC + 0.01),
              list(b0 = b0 - 0.05,
                   b = b - 0.5,
                   sig.quad = sig.quad + 0.04,
                   sig.transect = sig.transect + 0.04,
                   deltaA = deltaA + 0.04,
                   deltaB = deltaB + 0.04,
                   deltaC = deltaC + 0.04))

# JAGS model --------------------------------------------------------------

mod <- jagsUI::jags(data = data_list,
                    inits = inits,
                    #inits = NULL,
                    model.file = '/scratch/atm234/nps_plants/SAM/raw/inputs/plant_SAM_old_raw.R',
                    parameters.to.save = params,
                    parallel = TRUE,
                    n.chains = 3,
                    n.burnin = 5000,
                    n.iter =  45000,
                    n.thin = 10,
                    DIC = TRUE)

#save as an R data object
saveRDS(mod, 
        file ="/scratch/atm234/nps_plants/SAM/raw/outputs/nps_SAM_model2.RDS")

(end.time <- Sys.time())

(tot.time <- end.time - start.time)

# Check convergence -------------------------------------------------------

mcmcplot(mod$samples,
         dir = "/scratch/atm234/nps_plants/SAM/raw/outputs/mcmcplots/SAM2")

# Get RHat per parameter ------------------------------------------------

Rhat <- mod$Rhat

saveRDS(Rhat, "/scratch/atm234/nps_plants/SAM/raw/outputs/nps_SAM_model_Rhat2_raw.RDS")


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





