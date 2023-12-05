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
                     'fish_MSAM',
                     'outputs',
                     'fish_MSAM_model_RhatRE.RDS'))

rhat2 <- readRDS(here('01_sbc_fish',
                     'monsoon',
                     'fish_MSAM',
                     'outputs',
                     'fish_MSAM_modelRE_Rhat2.RDS'))


rhat4 <- readRDS(here('01_sbc_fish',
                      'monsoon',
                      'fish_MSAM',
                      'outputs',
                      'fish_MSAM_modelRE_Rhat4.RDS'))

rhat5 <- readRDS(here('01_sbc_fish',
                      'monsoon',
                      'fish_MSAM',
                      'outputs',
                      'fish_MSAM_modelRE_Rhat5.RDS'))

rhat6 <- readRDS(here('01_sbc_fish',
                      'monsoon',
                      'fish_MSAM',
                      'outputs',
                      'fish_MSAM_modelRE_Rhat6.RDS'))

# Graph RHat per parameter ------------------------------------------------

rhat_graph_fun(list = rhat, rhat = 1.2)

rhat_graph_fun(list = rhat2, rhat = 1.2)

rhat_graph_fun(rhat4)

rhat_graph_fun(rhat5)

parms <- c("a1.Vis", "a2.Size",  "eps.site.star",
           'eps.year.star', 'mu.a0', 'mu.b0species',
           'sig.a0', 'sig.b0species', 'sig.eps.site', 'sig.eps.year',
           "deviance", 'b0.star')

rhat_graph_fun(list = rhat6, parms = parms, rhat = 1.2) + #this is the one I'm using!!
labs(title = "SBC LTER fish MSAM Rhat")

ggsave(plot = last_plot(),
       filename = here("pictures",
                       "detection_models",
                       "SBC_Rhat_graph.jpg"),
       height = 7,
       width = 8,
       units = "in")

# summary of Rhat ---------------------------------------------------------

rhats <- rhat_df_fun(rhat6)

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

