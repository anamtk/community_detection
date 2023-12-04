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

source(here('00_functions',
            'plot_functions.R'))

# Load Data ---------------------------------------------------------------

rhat <- readRDS(here('02_konza_birds',
                     'monsoon',
                     "MSAM",
                     'outputs',
                     'bird_MSAM_model_RE_Rhat.RDS'))

rhat2 <- readRDS(here('02_konza_birds',
                     'monsoon',
                     "MSAM",
                     'outputs',
                     'bird_MSAM_modelRE_Rhat2.RDS'))


# Graph RHat per parameter ------------------------------------------------

parm <- c("a1.Effort", "a2.Size", 'b0.star', 'deviance',
          'eps.site.star', 'eps.year.star', 'mu.b0species',
          'sig.b0species', 'sig.eps.site', 'sig.eps.year',
          'mu.a0', 'sig.a0')

rhat_graph_fun(rhat, parms = parm, rhat = 1.2) + #this is the one I'm using!!
  labs(title = "KNZ LTER bird MSAM Rhat")

ggsave(plot = last_plot(),
       filename = here("pictures",
                       "detection_models",
                       "KNZ_Rhat_graph.jpg"),
       height = 7,
       width = 8,
       units = "in")

#rhat_graph_fun(rhat2)


# Rhat summaries ----------------------------------------------------------

rhats <- rhat_df_fun(rhat)

#remove ones we don't need to track
rhats2 <- rhats %>%
  filter(!id %in% c("eps.year", "eps.site"))

a <- rhats2 %>%
  filter(Rhat > 1.2) %>%
  tally()

b <- rhats2 %>%
  filter(Rhat <= 1.2) %>%
  tally()

1-(a/(a+b))

rhats_root <- rhats2 %>%
  filter(!id == "b0.star")

c <- rhats_root %>%
  filter(Rhat > 1.2) %>%
  tally()

d <- rhats_root %>%
  filter(Rhat <= 1.2) %>%
  tally()

c/(c+d)
