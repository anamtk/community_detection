#Getting initials for next model run
#Ana Miller-ter Kuile
#September 15, 2023

#this script pulls out the last values from the chain with the lowest deviance
#to use as initials in a follow-up model run



# Load packages ---------------------------------------------------------------
Sys.time()


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

mod <- readRDS( file ="/scratch/atm234/sev_hoppers/outputs/sev_MSAM_nocov_model.RDS")

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

#for bird model root nodes:
#omega
mu.llambda <- as.vector(samp_df2$mu.llambda)
sig.llambda <- as.vector(samp_df2$sig.llambda)
p.mean <- as.vector(samp_df2$p.mean)
sig.lp <- as.vector(samp_df2$sig.lp)


# Load Data ---------------------------------------------------------------

#load the formatted data for the JAGS model
data <- readRDS("/scratch/atm234/sev_hoppers/inputs/sev_msam_dynmultisite.RDS")

# Compile data ------------------------------------------------------------

data_list <-  list(y = data$y,
                   reprod = data$reprod,
                   n.species = data$n.species,
                   n.years = data$n.years,
                   n.start = data$n.start,
                   n.end = data$n.end,
                   n.transects = data$n.transects,
                   n.rep = data$n.rep,
                   #for initials
                   ymax = data$ymax,
                   omega.init = data$omega.init,
                   #for omega prior
                   R = data$R)



# Parameters to save ------------------------------------------------------


params <- c(
  #COMMUNITY parameters
  'mu.llambda',
  'sig.llambda',
  'p.mean',
  'sig.lp')

# INits -------------------------------------------------------------------

#we found ymax to set initials, since otherwise the model will hate us
#also Kiona suggested setting initials for omega based on covariance, since
#the model will struggle with this
inits <- list(list(N = data$ymax,
                              mu.llambda = mu.llambda,
                              sig.llambda = sig.llambda,
                              p.mean = p.mean,
                              sig.lp = sig.lp),
                         list(N = data$ymax,
                              mu.llambda = mu.llambda + 0.05,
                              sig.llambda = sig.llambda + 0.04,
                              p.mean = p.mean + 0.005,
                              sig.lp = sig.lp + 0.05),
                         list(N = data$ymax,
                              mu.llambda = mu.llambda - 0.05,
                              sig.llambda = sig.llambda + 0.08,
                              p.mean = p.mean -0.005,
                              sig.lp = sig.lp + 0.08)
                         )

# JAGS model --------------------------------------------------------------

mod2 <- jagsUI::jags(data = data_list,
                    inits = inits,
                    #inits = NULL,
                    model.file = '/scratch/atm234/sev_hoppers/inputs/sev_MSAM_simple_nocov.R',
                    parameters.to.save = params,
                    parallel = TRUE,
                    n.chains = 3,
                    n.iter = 100000,
                    n.burnin = 10000,
                    n.thin = 20,
                    DIC = TRUE)

#save as an R data object
saveRDS(mod2,
        file ="/scratch/atm234/sev_hoppers/outputs/sev_MSAMnocov_model2.RDS")

Sys.time()



# Check convergence -------------------------------------------------------

mcmcplot(mod2$samples,
         dir = "/scratch/atm234/sev_hoppers/outputs/mcmcplots/MSAM2_nocov")

# Get RHat per parameter ------------------------------------------------

Rhat <- mod2$Rhat

saveRDS(Rhat, "/scratch/atm234/sev_hoppers/outputs/sev_MSAMnocov_model_Rhat2.RDS")


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





