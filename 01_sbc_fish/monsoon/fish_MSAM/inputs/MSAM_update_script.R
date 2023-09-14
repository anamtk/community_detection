#fish MSAM update
#Ana Miller-ter Kuile
#September 7, 2023

#this script updates the MSAM based on raftery of initial run


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

# Load Model ---------------------------------------------------------------

#load the model
mod <- readRDS(file ="/scratch/atm234/sbc_fish/outputs/fish_MSAM_model.RDS")



# Parameters to save ------------------------------------------------------


params <- c(
  #COMMUNITY parameters
  'a1.Vis',
  'a2.Size',
  'lambda.mean',
  'sig.lambda',
  'a0.mean',
  'sig.a0',
  'omega')

# Update model ------------------------------------------------------------
#31590

mod2 <- update(mod,
               parameters.to.save = params,
               n.iter = 31590,
               parallel = T)

#save as an R data object
saveRDS(mod2, 
        file ="/scratch/atm234/sbc_fish/outputs/fish_MSAM_model_update.RDS")

Sys.time()



# Check convergence -------------------------------------------------------

mcmcplot(mod2$samples,
         dir = "/scratch/atm234/sbc_fish/outputs/mcmcplots/MSAM/update")

# Get RHat per parameter ------------------------------------------------

Rhat <- mod2$Rhat

saveRDS(Rhat, "/scratch/atm234/sbc_fish/outputs/fish_MSAM_model_Rhat_update.RDS")


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




