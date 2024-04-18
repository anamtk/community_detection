#Convergence check for nps SAM
#Shelby Lamm
#November 30, 2023

#this script assesses convergence in the nps sam

# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 
package.list <- c("here", "tidyverse") 


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

source(here('00_functions',
            'plot_functions.R'))

# Load Data ---------------------------------------------------------------

rhat <- readRDS(here('04_nps_plants',
                     'monsoon',
                     "nps_SAM",
                     'outputs',
                     'nps_SAM_model_Rhat.RDS'))


# Graph RHat per parameter ------------------------------------------------

#rhat_graph_fun(rhat)

parm <- c('b0.quad',
            'b0.transect',
            'b',
            'wA',
            'wB',
            'sig.quad',
            'sig.transect',
            'var.process', 
          'deviance')

rhat_graph_fun(rhat, parms = parm, rhat = 1.1) +
  labs(title = "NPS Plant SAM Rhat")

ggsave(plot = last_plot(),
       filename = here("pictures",
                       "supplementary",
                       'SAM',
                       "NPS_SAM_Rhat_graph.jpg"),
       height = 4,
       width = 6,
       units = "in")


