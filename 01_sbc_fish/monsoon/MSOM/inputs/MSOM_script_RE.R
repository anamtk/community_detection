#Monsoon script - fish MSAM
# Ana Miller-ter Kuile
# July 27, 2023

#this script runs the fish MSOM

# Load packages ---------------------------------------------------------------
(start.time <- Sys.time())


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
data <- readRDS("/scratch/atm234/sbc_fish/inputs/fish_msom_RE.RDS")

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
                   n.sites = data$n.sites,
                   Site.ID = data$Site.ID,
                   Year.ID = data$Year.ID,
                   #for initials
                   z = data$z)

# Parameters to save ------------------------------------------------------

params <- c(
  #COMMUNITY parameters
  'b0.star',
  'eps.site.star',
  'eps.year.star',
  'eps.site',
  'eps.year',
  'a1.Vis',
  'a2.Size',
  'mu.a0',
  'sig.a0',
  'mu.b0species',
  'sig.b0species',
  'sig.eps.site',
  'sig.eps.year'
)

# INits -------------------------------------------------------------------

#we found ymax to set initials, since otherwise the model will hate us
# inits <- list(list(N = data$ymax),
#               list(N = data$ymax),
#               list(N = data$ymax))

# JAGS model --------------------------------------------------------------

mod <- jagsUI::jags(data = data_list,
                    #inits = inits,
                    inits = NULL,
                    model.file = '/scratch/atm234/sbc_fish/inputs/MSOM_simple_siteyearRE.R',
                    parameters.to.save = params,
                    parallel = TRUE,
                    n.chains = 3,
                    #n.burnin = 1000,
                    #n.thin = 3,
                    n.iter = 4000,
                    DIC = TRUE)

#save as an R data object
saveRDS(mod, 
        file ="/scratch/atm234/sbc_fish/outputs/fish_MSOM_model_RE.RDS")

(end.time <- Sys.time())


(tot.time <- end.time - start.time)
# Check convergence -------------------------------------------------------

params2 <- c(
  #COMMUNITY parameters
  'b0.star',
  'eps.site.star',
  'eps.year.star',
  'a1.Vis',
  'a2.Size',
  'mu.a0',
  'sig.a0',
  'mu.b0species',
  'sig.b0species',
  'sig.eps.site',
  'sig.eps.year',
  'deviance'
)

mcmcplot(mod$samples,
         parms = params2,
         dir = "/scratch/atm234/sbc_fish/outputs/mcmcplots/MSOM_RE")

# Get RHat per parameter ------------------------------------------------

Rhat <- mod$Rhat

saveRDS(Rhat, "/scratch/atm234/sbc_fish/outputs/fish_MSOM_model_RhatRE.RDS")


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




