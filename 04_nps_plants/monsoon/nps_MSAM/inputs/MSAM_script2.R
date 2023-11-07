#Getting initials for next model run
#Shelby Lamm
#October 24, 2023

#this script pulls out the last values from the chain with the lowest deviance
#to use as initials in a follow-up model run

# Load packages ---------------------------------------------------------------
(start.time <- Sys.time())


# Load packages,
package.list <- c("jagsUI", "coda",
                  'dplyr', 'stringr',
                  'magrittr', 'tidyr',
                  'mcmcplots','ggplot2',
                  'tibble') 


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load model --------------------------------------------------------------

mod <- readRDS(file ="/scratch/sml665/nps_plants/outputs/nps_MSAM_model.RDS")

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


# for nps model root nodes:
mu.lpsi <- as.vector(samp_df2$mu.lpsi)
sig.lpsi <- as.vector(samp_df2$sig.lpsi)
mu.a0 <- as.vector(samp_df2$mu.a0)
sig.a0 <- as.vector(samp_df2$sig.a0)
a1.Cover <- as.vector(samp_df2$a1.Cover)
a2.Lifegroup <- as.vector(samp_df2$a2.Lifegroup)
#mu.missingcover = as.vector(samp_df2$mu.missingcover)
#sig.missingcover = as.vector(samp_df2$sig.missingcover)


# Load Data ---------------------------------------------------------------

#load the formatted data for the JAGS model
#data <- readRDS("/scratch/sml665/nps_plants/inputs/nps_msam_multisite.RDS")
data <- readRDS("/scratch/sml665/nps_plants/inputs/nps_msam_multisite_subset.RDS")

# Compile data ------------------------------------------------------------
data_list <- list(n.species = data$n.species,
                  n.quads = data$n.quads,
                  n.yr = data$n.yr,
                  n.rep = data$n.rep,
                  cover = data$cover,
                  lifegroup = data$lifegroup,
                  n.groups = data$n.groups,
                  y = data$y,
                  z = data$z)

# Parameters to save ------------------------------------------------------

params <- c(
  #COMMUNITY parameters
  'mu.lpsi',
  'sig.lpsi',
  'mu.a0',
  'sig.a0',
  'a1.Cover',
  'a2.LifeGroup' #,
  #'p0',
  #'mu.missingcover',
  #'sig.missingcover'
)


# INits -------------------------------------------------------------------

#we found z to set initials, since otherwise the model will hate us
inits <- list(list(N = data$z,
                   mu.lpsi = samp_df2$mu.lpsi,
                   sig.lpsi = samp_df2$sig.lpsi,
                   mu.a0 = samp_df2$mu.a0,
                   sig.a0 = samp_df2$sig.a0,
                   a1.Cover = samp_df2$a1.Cover,
                   a2.Lifegroup = a2.Lifegroup),
                   #mu.missingcover = samp_df2$mu.missingcover,
                   #sig.missingcover = samp_df2$sig.missingcover),
              list(N = data$z,
                   mu.lpsi = samp_df2$mu.lpsi + 0.25,
                   sig.lpsi = samp_df2$sig.lpsi + 0.05,
                   mu.a0 = samp_df2$mu.a0 + 0.25,
                   sig.a0 = samp_df2$sig.a0 + 0.05,
                   a1.Cover = samp_df2$a1.Cover + 0.02,
                   a2.Lifegroup = a2.Lifegroup + 0.4),
                   #mu.missingcover = samp_df2$mu.missingcover + 0.003,
                   #sig.missingcover = samp_df2$sig.missingcover + 0.003),
              list(N = data$z,
                   mu.lpsi = samp_df2$mu.lpsi - 0.25,
                   sig.lpsi = samp_df2$sig.lpsi + 0.1,
                   mu.a0 = samp_df2$mu.a0 - 0.25,
                   sig.a0 = samp_df2$sig.a0 + 0.1,
                   a1.Cover = samp_df2$a1.Cover - 0.02,
                   a2.Lifegroup = a2.Lifegroup - 0.4)) #,
                   #mu.missingcover = samp_df2$mu.missingcover - 0.003,
                   #sig.missingcover = samp_df2$sig.missingcover - 0.003))

# JAGS model --------------------------------------------------------------

mod2 <- jagsUI::jags(data = data_list,
                    inits = inits,
                    #inits = NULL,
                    model.file = '/scratch/sml665/nps_plants/inputs/nps_MSOM_simple_AMtK.R',
                    parameters.to.save = params,
                    parallel = TRUE,
                    n.chains = 3,
                    n.iter = 50000,
                    n.burnin = 5000,
                    n.thin = 10,
                    DIC = TRUE)

#save as an R data object
saveRDS(mod2, 
        file ="/scratch/sml665/nps_plants/outputs/nps_MSAM_model2.RDS")

(end.time <- Sys.time())


(tot.time <- end.time - start.time)
# Check convergence -------------------------------------------------------

mcmcplot(mod2$samples,
         dir = "/scratch/sml665/nps_plants/outputs/mcmcplots/MSAM2")

# Get RHat per parameter ------------------------------------------------

Rhat <- mod2$Rhat

saveRDS(Rhat, "/scratch/sml665/nps_plants/outputs/nps_MSAM_model_Rhat2.RDS")


# Get Raftery diag --------------------------------------------------------


raf <- raftery.diag(mod2$samples)

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






