# Running the MSAM for fish
# Ana Miller-ter Kuile
# July 25, 2023

# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 
package.list <- c("here", "tidyverse", 
                  "jagsUI",
                  'rjags',
                  'mcmcplots',
                  "coda") #mcmc output


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load Data ---------------------------------------------------------------

data <- readRDS(here('01_sbc_fish',
                     "data_outputs",
                     'MSAM',
                     "model_inputs",
                     "fish_msam_dynmultisiteRE.RDS"))

data <- readRDS(here('01_sbc_fish',
                     "data_outputs",
                     'MSAM',
                     "model_inputs",
                     "fish_msom_RE.RDS"))

# Parameters to save ------------------------------------------------------

params <- c(
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
            'eps.site',
            'eps.year'
            )

#params <- c("bray")

#we found ymax to set initials, since otherwise the model will hate us
# inits <- list(list(N = data$ymax,
#                    omega = data$omega.init),
#               list(N = data$ymax,
#                    omega = data$omega.init),
#               list(N = data$ymax,
#                    omega = data$omega.init))

inits <- list(list(N = data$ymax),
              list(N = data$ymax),
              list(N = data$ymax))

inits <- list(list(z = data$z),
              list(z = data$z),
              list(z = data$z))

#inits <- function()list(N = data$ymax)

# JAGS model --------------------------------------------------------------

model <- here("01_sbc_fish",
              "code", 
              "03_sb_analyses",
              '01_MSAM_modeling',
              "models",
              "MSOM_simple_siteyearRE.R")

(st.time <- Sys.time())
mod <- jagsUI::jags(data = data,
                    inits = NULL,
                    #inits = inits,
                    model.file = model,
                    parameters.to.save = params,
                    parallel = TRUE,
                    n.chains = 3,
                    #n.burnin = 2000,
                    n.iter = 1,
                    DIC = TRUE)

end.time <- Sys.time()

end.time - st.time
# 
#mcmcplot(mod$samples)
# 
#gelman.diag(mod$samples)
# 
# # Get inits ---------------------------------------------------------------
# 
# #get the MCMC chains
# samples <- mod$samples
# 
# #function to make each chain a dataframe
# df_fun <- function(chain){
#   df <- as.data.frame(chain) %>%
#     rownames_to_column(var = "iteration")
#   return(df)
# }
# 
# #use that function on all list elements
# samp_dfs <- lapply(samples, df_fun)
# 
# #make into one dataframe
# samp_df <- bind_rows(samp_dfs, .id = "chain")
# 
# #get values for all parameters from the last iteration of the
# #chain with the lowest deviance
# samp_df2 <- samp_df %>%
#   group_by(chain) %>%
#   #get mean deviance by chain
#   mutate(mean_dev = mean(deviance, na.rm = T)) %>%
#   ungroup() %>%
#   #get only the chain with the minimum average deviance
#   filter(mean_dev == min(mean_dev)) %>%
#   #pull out the final iteration from that chain
#   filter(iteration == max(iteration)) %>%
#   dplyr::select(-chain, -iteration,
#                 -deviance, -mean_dev) 
# 
# #for fish model root nodes:
# eps.site<- samp_df2 %>%
#   dplyr::select(contains('eps.site')) %>%
#   dplyr::select(-contains('eps.site.star')) %>%
#   dplyr::select(-contains('sig')) %>%
#   pivot_longer(cols = everything(),
#                names_to = "parm",
#                values_to = "value") %>%
#   separate(parm, 
#            into = c("site", "species"),
#            sep = ",") %>%
#   mutate(site = str_sub(site, 10, nchar(site))) %>%
#   mutate(species = str_sub(species, 1, (nchar(species)-1))) %>%
#   pivot_wider(names_from = species,
#               values_from = value) %>%
#   dplyr::select(-site) %>%
#   as.matrix()
# 
# eps.year<- samp_df2 %>%
#   dplyr::select(contains('eps.year')) %>%
#   dplyr::select(-contains('eps.year.star')) %>%
#   dplyr::select(-contains('sig')) %>%
#   pivot_longer(cols = everything(),
#                names_to = "parm",
#                values_to = "value") %>%
#   separate(parm, 
#            into = c("year", "species"),
#            sep = ",") %>%
#   mutate(year = str_sub(year, 10, nchar(year))) %>%
#   mutate(species = str_sub(species, 1, (nchar(species)-1))) %>%
#   pivot_wider(names_from = species,
#               values_from = value) %>%
#   dplyr::select(-year) %>%
#   as.matrix()
# 
# mu.a0 <- as.vector(samp_df2$mu.a0)
# sig.a0<- as.vector(samp_df2$sig.a0)
# a1.Vis <- as.vector(samp_df2$a1.Vis)
# a2.Size <- as.vector(samp_df2$a2.Size)
# mu.b0species <- as.vector(samp_df2$mu.b0species)
# sig.b0species <- as.vector(samp_df2$sig.b0species)
# sig.eps.site <- as.vector(samp_df2$sig.eps.site)
# sig.eps.year <- as.vector(samp_df2$sig.eps.year)
