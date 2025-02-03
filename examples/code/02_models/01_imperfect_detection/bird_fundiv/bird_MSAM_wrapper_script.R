# Running the MSAM for birds
# September 8, 2023

# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 
package.list <- c("here", "tidyverse", 
                  'data.table',
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

#for plotting convergence checks
source(here('examples',
            'code',
            '00_functions',
            'plot_functions.R'))


# Load Data ---------------------------------------------------------------

data <- readRDS(here('examples',
                     'data_output',
                     'bird_fundiv',
                     'model_inputs',
                     "bird_msam_data_input_list.RDS"))

# Parameters to save ------------------------------------------------------

params <- c(
            #COMMUNITY parameters
            'a1.Effort',
            'a2.Size',
            'b0.star',
            'eps.site.star',
            'eps.year.star',
            'mu.b0species',
            'sig.b0species',
            'sig.eps.site',
            'sig.eps.year',
            'mu.a0',
            'sig.a0')


#we found ymax to set initials, since otherwise the model will hate us
inits <- function() list(N = data$ymax)#,
                         #omega = data$omega.init)

# JAGS model --------------------------------------------------------------

model <- here('examples',
              'code',
              '02_models',
              '01_imperfect_detection',
              'bird_fundiv',
              'bird_MSAM_model.R')

#this model takes a long time to run and we parallelized it on a 
#computing cluster, where it took ~hours/days to run. 
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


# Convergence and GOF -----------------------------------------------------

rhat <- mod$Rhat
#once the model runs, you can evaluate convergence and calculate GOF 
#for a converged model
parm <- c("a1.Effort", "a2.Size", 'b0.star', 'deviance',
          'eps.site.star', 'eps.year.star', 'mu.b0species',
          'sig.b0species', 'sig.eps.site', 'sig.eps.year',
          'mu.a0', 'sig.a0')

rhat_graph_fun(rhat, parms = parm, rhat = 1.2) + #this is the one I'm using!!
  labs(title = "KNZ LTER bird MSAM Rhat")

# GOF ---------------------------------------------------------------------


# Load data ---------------------------------------------------------------

observed <- read.csv(here('examples',
                          "data_output",
                          'bird_fundiv',
                          'tidy_data',
                          'bird_msam_tidy_data.csv'))

#again, generated from the computing cluster after
#model convergence
modeled <- readRDS(here('examples',
                        "data_output",
                        'bird_fundiv',
                        'computing_cluster_outputs',
                        'bird_MSAM_GOF_summary.RDS'))


# Pull out yrep -----------------------------------------------------------

#y.rep is indexed species, transects, years, replicates

yrep <- as.data.frame(modeled$statistics) %>%
  rownames_to_column(var = 'parm') %>%
  filter(str_detect(parm, "y.rep")) %>%
  separate(parm,
           into = c("SpecID", "TransID", "yrID", "REP"),
           sep = ",") %>%
  mutate(SpecID = str_sub(SpecID, 7, nchar(SpecID))) %>%
  mutate(REP = str_sub(REP, 1, (nchar(REP)-1))) %>%
  mutate(SpecID = as.numeric(SpecID),
         TransID = as.numeric(TransID),
         yrID = as.numeric(yrID),
         REP = as.numeric(REP)) %>%
  left_join(observed, by = c("SpecID", "TransID", "yrID", "REP"))


# Plot --------------------------------------------------------------------

m1 <- lm(Mean ~ NOBS,
         data = yrep)

summary(m1)

lb1 <- paste("R^2 == 0.48")

ggplot(yrep, aes(x = NOBS, y = Mean)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point()+
  annotate(geom = 'text', x = 125, y = 20, label = lb1, parse = T) +
  labs(x = "observed", y = "predicted", title = "KNZ LTER bird MSAM GOF")

