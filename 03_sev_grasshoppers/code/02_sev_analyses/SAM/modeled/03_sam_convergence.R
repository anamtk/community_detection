#Convergence check for fish MSAM
# Ana Miller-ter Kuile
#Sepptember 7, 2023

#this script assesses convergence in the fish msam

# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 
package.list <- c("here", "tidyverse") 


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

source(here("00_functions",
            'plot_functions.R'))

# Load Data ---------------------------------------------------------------

rhat <- readRDS(here('03_sev_grasshoppers',
                     'monsoon',
                     'SAM',
                     'outputs',
                     'sev_SAM_model_Rhat.RDS'))

rhat2 <- readRDS(here('03_sev_grasshoppers',
                     'monsoon',
                     'SAM',
                     'outputs',
                     'sev_SAM_model_Rhat2.RDS'))


rhat4 <- readRDS(here('03_sev_grasshoppers',
                      'monsoon',
                      'SAM',
                      'outputs',
                      'sev_SAM_model_Rhat4.RDS'))

# Graph RHat per parameter ------------------------------------------------

parms <- c('b0.web', 'b0.transect',
           'b', 'b0', 'wA', 'wB', "wC",
           'sig.web', 'sig.transect', 'var.process',
           'deviance')

rhat_graph_fun(rhat2, parms = parms, rhat = 1.1)

rhat_graph_fun(rhat4, parms = parms, rhat = 1.1 ) +
  labs(title = "SEV LTER fish SAM Rhat")

ggsave(plot = last_plot(),
       filename = here("pictures",
                       "supplementary",
                       'SAM',
                       "SEV_SAM_Rhat_graph.jpg"),
       height = 4,
       width = 6,
       units = "in")
