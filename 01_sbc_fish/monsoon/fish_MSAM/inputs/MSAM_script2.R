#Getting initials for next model run
#Ana Miller-ter Kuile
#September 15, 2023

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

mod <- readRDS(file ="/scratch/atm234/sbc_fish/outputs/fish_MSAM_model.RDS")

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

#for fish model root nodes:
mu.llambda <- as.vector(samp_df2$mu.llambda)
sig.llambda <- as.vector(samp_df2$sig.llambda)
mu.a0 <- as.vector(samp_df2$mu.a0)
#sig.a0 seemed to be covering the whole domain of dunif(0,50) so updated to 100 in
#next run
sig.a0 <- as.vector(samp_df2$sig.a0)
a1.Vis <- as.vector(samp_df2$a1.Vis)
a2.Size <- as.vector(samp_df2$a2.Size)


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
  'mu.llambda',
  'sig.llambda',
  'mu.a0',
  'sig.a0'
)


# INits -------------------------------------------------------------------

#we found ymax to set initials, since otherwise the model will hate us
inits <- list(list(N = data$ymax,
                   mu.llambda = mu.llambda,
                   sig.llambda = sig.llambda,
                   mu.a0 = mu.a0,
                   sig.a0 = sig.a0,
                   a1.Vis = a1.Vis,
                   a2.Size = a2.Size),
              list(N = data$ymax,
                   mu.llambda = mu.llambda + 0.05,
                   sig.llambda = sig.llambda + 0.05,
                   mu.a0 = mu.a0 + 0.05,
                   sig.a0 = sig.a0 + 0.05,
                   a1.Vis = a1.Vis + 0.05,
                   a2.Size = a2.Size + 0.5),
              list(N = data$ymax,
                   mu.llambda = mu.llambda - 0.05,
                   sig.llambda = sig.llambda + 1,
                   mu.a0 = mu.a0 - 0.05,
                   sig.a0 = sig.a0 + 0.08,
                   a1.Vis = a1.Vis - 0.05,
                   a2.Size = a2.Size - 0.5))

# JAGS model --------------------------------------------------------------

mod2 <- jagsUI::jags(data = data_list,
                    inits = inits,
                    #inits = NULL,
                    model.file = '/scratch/atm234/sbc_fish/inputs/MSAM_simple.R',
                    parameters.to.save = params,
                    parallel = TRUE,
                    n.chains = 3,
                    n.iter = 100000,
                    n.burnin = 10000,
                    n.thin = 20,
                    DIC = TRUE)

#save as an R data object
saveRDS(mod2, 
        file ="/scratch/atm234/sbc_fish/outputs/fish_MSAM_model2.RDS")

(end.time <- Sys.time())


(tot.time <- end.time - start.time)
# Check convergence -------------------------------------------------------

mcmcplot(mod2$samples,
         dir = "/scratch/atm234/sbc_fish/outputs/mcmcplots/MSAM2")

# Get RHat per parameter ------------------------------------------------

Rhat <- mod2$Rhat

saveRDS(Rhat, "/scratch/atm234/sbc_fish/outputs/fish_MSAM_model_Rhat2.RDS")


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






