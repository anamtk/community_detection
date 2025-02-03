# Running the MSAM for grasshoppers
# September 11, 2023

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

source(here('examples',
            'code',
            '00_functions',
            'plot_functions.R'))

# Load Data ---------------------------------------------------------------

data <- readRDS(here('examples',
                     'data_output',
                     'grasshopper_stability',
                     'model_inputs',
                     "grasshopper_msam_input_data_list.RDS"))

# Parameters to save ------------------------------------------------------

params <- c(
            #COMMUNITY parameters
            'b0.star',
            'eps.site.star',
            'eps.year.star',
            'p.mean',
            'sig.lp', 
            'mu.b0species',
            'sig.b0species',
            'sig.eps.site',
            'sig.eps.year')


#we found ymax to set initials, since otherwise the model will hate us
inits <- function() list(N = data$ymax)

# JAGS model --------------------------------------------------------------

model <- here('examples',
              'code',
              '02_models',
              '01_imperfect_detection',
              'grasshopper_stability',
              "grasshopper_MSAM_model.RDS")

#note that this model takes a very long time to run and 
#converge - so we ran these models on a computer cluster, often
#re-running the model with initial values for specific parameters
#to get a converged model to work with.

start.time <- Sys.time()
mod <- jagsUI::jags(data = data,
                         inits = inits,
                         model.file = model,
                         parameters.to.save = params,
                         parallel = TRUE,
                         n.chains = 3,
                         n.iter = 1,
                         DIC = TRUE)

end.time <- Sys.time()

end.time - start.time


# Check convergence -------------------------------------------------------

mcmcplot(mod$samples)
gelman.diag(mod$samples, multivariate = F)

rhat <- mod$Rhat

parm <- c("b0.star", 'eps.site.star', 'eps.year.star', 'p.mean',
          'sig.lp', 'mu.b0species','sig.b0species',
          'sig.eps.year', 'sig.eps.site')
rhat_graph_fun(list = rhat, parms = parm, rhat = 1.2)  + #this is the one I'm using!!
  labs(title = "SEV LTER grasshopper MSAM Rhat")

# GOF ---------------------------------------------------------------------

observed <- read.csv(here('examples',
                          'data_output',
                          'grasshopper_stability',
                          'tidy_data',
                          'grasshopper_msam_tidy_data.csv'))


modeled <- readRDS(here('examples',
                        'data_output',
                        'grasshopper_stability',
                        'computing_cluster_outputs',
                        'grasshopper_MSAM_GOF_summary.RDS'))

# Pull out yrep -----------------------------------------------------------

#y.rep is indexed species, transects, years, replicates

yrep <- as.data.frame(modeled$statistics) %>%
  rownames_to_column(var = 'parm') %>%
  filter(str_detect(parm, "y.rep")) %>%
  separate(parm,
           into = c("speciesID", "siteID", "yrID", "rep"),
           sep = ",") %>%
  mutate(speciesID = str_sub(speciesID, 7, nchar(speciesID))) %>%
  mutate(rep = str_sub(rep, 1, (nchar(rep)-1))) %>%
  mutate(speciesID = as.numeric(speciesID),
         siteID = as.numeric(siteID),
         yrID = as.numeric(yrID),
         rep = as.numeric(rep)) %>%
  left_join(observed, by = c("speciesID", "siteID", "yrID", "rep"))


# Plot --------------------------------------------------------------------

m1 <- lm(Mean ~ CNT,
         data = yrep)

summary(m1)

lb1 <- paste("R^2 == 0.45")

ggplot(yrep, aes(x = CNT, y = Mean)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point() +
  annotate(geom = 'text', x = 20, y = 2.5, label = lb1, parse = T) +
  labs(x = "observed", y = "predicted", title = "SEV LTER grasshopper MSAM GOF")





