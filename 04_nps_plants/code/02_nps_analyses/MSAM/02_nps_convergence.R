#Convergence check for nps MSAM
#Shelby Lamm
#October 30, 2023

#this script assesses convergence in the nps msam

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
                     "nps_MSAM",
                     'outputs_yrsite',
                     'nps_MSAM_model_RE_Rhat.RDS'))


# Graph RHat per parameter ------------------------------------------------

parm <- c("b0.star", "eps.site.star", "eps.year.star",
          'sig.eps.site', 'sig.eps.year', 'deviance',
          'mu.b0species', 'sig.b0species',
          'mu.a0', 'sig.a0','a1.Cover', 'a2.LifeGroup')

rhat_graph_fun(rhat, parms = parm, rhat = 1.2)+ 
  labs(title = "PFNP plant MSOM Rhat")

ggsave(plot = last_plot(),
       filename = here("pictures",
                       "detection_models",
                       "PFNP_Rhat_graph.jpg"),
       height = 7,
       width = 8,
       units = "in")

