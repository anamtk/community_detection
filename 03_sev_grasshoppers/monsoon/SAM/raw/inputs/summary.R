#Get Bray-Curtis out of converged model
#Ana Miller-ter Kuile
#September 19, 2023

#this script pulls out and sumamrises Bray Curtis for visualizations
#purposes for An's project


# Load packages ---------------------------------------------------------------
Sys.time()


# Load packages,
package.list <- c('jagsUI',"coda", 'dplyr',
                  'stringr', 'purrr', 'tibble',
                  'magrittr', 'tidyr') 


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load model --------------------------------------------------------------

mod <- readRDS( file = "/scratch/atm234/sev_hoppers/SAM/raw/outputs/sev_SAM_model_raw.RDS")


# Get summary of model -------------------------------------------------

sum <- summary(mod$samples)

saveRDS(sum, 
        file = "/scratch/atm234/sev_hoppers/SAM/raw/outputs/sev_SAM_raw_summary.RDS")


# cummulative weights -----------------------------------------------------

samples <- mod$sims.list

wA <- as.data.frame(samples$wA) %>%
  rownames_to_column(var = 'iteration') %>%
  pivot_longer(V1:V6,
               names_to = 'weightID',
               values_to = "value") %>%
  group_by(iteration) %>%
  mutate(cumm.weight = accumulate(value, `+`)) %>%
  ungroup() %>%
  dplyr::select(weightID, cumm.weight) %>%
  group_by(weightID) %>%
  summarise(median_obs = median(cumm.weight),
            LCI_obs = quantile(cumm.weight, 0.025),
            UCI_obs = quantile(cumm.weight, 0.975)) %>%
  mutate(lag = substr(weightID, 2, (nchar(weightID)))) %>%
  mutate(lag = as.numeric(lag)) %>%
  mutate(variable = "Temperature") %>%
  mutate(dataset = "SEV grasshoppers") %>%
  dplyr::select(lag, variable, dataset,
                LCI_obs, median_obs, UCI_obs)

saveRDS(wA, 
        file = "/scratch/atm234/sev_hoppers/SAM/raw/outputs/sev_SAM_cummtemp_raw.RDS")

wB <- as.data.frame(samples$wB) %>%
  rownames_to_column(var = 'iteration') %>%
  pivot_longer(V1:V6,
               names_to = 'weightID',
               values_to = "value") %>%
  group_by(iteration) %>%
  mutate(cumm.weight = accumulate(value, `+`)) %>%
  ungroup() %>%
  dplyr::select(weightID, cumm.weight) %>%
  group_by(weightID) %>%
  summarise(median_obs = median(cumm.weight),
            LCI_obs = quantile(cumm.weight, 0.025),
            UCI_obs = quantile(cumm.weight, 0.975))%>%
  mutate(lag = substr(weightID, 2, (nchar(weightID)))) %>%
  mutate(lag = as.numeric(lag)) %>%
  mutate(variable = "Precipitation") %>%
  mutate(dataset = "SEV grasshoppers") %>%
  dplyr::select(lag, variable, dataset,
                LCI_obs, median_obs, UCI_obs)

saveRDS(wB, 
        file = "/scratch/atm234/sev_hoppers/SAM/raw/outputs/sev_SAM_cummppt_raw.RDS")
