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

rhat <- readRDS(here('01_sbc_fish',
                     'monsoon',
                     'SAM',
                     'outputs',
                     'fish_SAM_model_Rhat.RDS'))



# Graph RHat per parameter ------------------------------------------------
parm <- c('b0.transect', 'b', 'wA', 'wB', 'sig.transect', 'sig.site',
          'var.process', 'deviance')

rhat_graph_fun(rhat, parms = parm, rhat = 1.1) +
  labs(title = "SBC LTER fish SAM Rhat")

ggsave(plot = last_plot(),
       filename = here("pictures",
                       "supplementary",
                       'SAM',
                       "SBC_SAM_Rhat_graph.jpg"),
       height = 4,
       width = 6,
       units = "in")
