#Monsoon script - fish MSAM
# Ana Miller-ter Kuile
# July 27, 2023

#this script runs the fish MSOM

# Load packages ---------------------------------------------------------------
Sys.time()


# Load packages,
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
data <- readRDS("nps_msam_dynmultisite.RDS")

# Compile data ------------------------------------------------------------

data_list <- list(n.species = data$n.species,
             n.quads = data$n.quads,
             n.yr = data$n.yr,
             n.rep = data$n.rep,
             y = data$y,
             #omega prior:
             R = data$R,
             #initials
             z = data$z,
             omega.init1 = data$omega.init1,
             omega.init2 = data$omega.init2,
             omega.init3 = data$omega.init3)

# Parameters to save ------------------------------------------------------

params <- c(
  #COMMUNITY parameters
  'psi.mean',
  'sig.lpsi',
  'p.mean',
  'sig.lp',
  'omega'
)

# INits -------------------------------------------------------------------

#we found ymax to set initials, since otherwise the model will hate us
#also Kiona suggested setting initials for omega based on covariance, since
#the model will struggle with this

#inits <- function() list(N = data$z) #,

#model breaks with these initials...
                         #omega = data$omega.init)

inits <- list(list(N = data$z,
                   omega = data$omega.init1),
              list(N = data$z,
                   omega = data$omega.init2),
              list(N = data$z,
                   omega = data$omega.init3))

# JAGS model --------------------------------------------------------------

mod <- jagsUI::jags(data = data_list,
                        inits = inits,
                        #inits = NULL,
                        model.file = 'nps_dyn_MSOM_cov.R',
                        parameters.to.save = params,
                        parallel = TRUE,
                        n.chains = 3,
                        n.iter = 4000,
                        DIC = TRUE)

#save as an R data object
saveRDS(mod, 
        file ="./outputs/nps_MSAM_model.RDS")

Sys.time()



# Check convergence -------------------------------------------------------

#i don't make MCMC plots of omega because it's a HUGE parameter
parms <- c(
  #COMMUNITY parameters
  'psi.mean',
  'sig.lpsi',
  'p.mean',
  'sig.lp'
)

mcmcplot(mod$samples,
         parms = parms,
         dir = "./outputs/mcmcplots/MSAM")

# Get RHat per parameter ------------------------------------------------

#I do pull omega out here - and use it in an R script on my computer
# to plot per-parameter convergence
Rhat <- mod$Rhat

saveRDS(Rhat, "./outputs/nps_MSAM_model_Rhat.RDS")


# Get Raftery diag --------------------------------------------------------

#if you look at your .out file on monsoon - you can
#see what Raftery thinks you need to run in terms of 
#iterations and burnin from the following code
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



