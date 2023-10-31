#Detection partial effects plots
#Ana Miller-ter Kuile
#October 27, 2023

#this script creates partial plots of the detection covariates in the MSAMs

# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 
package.list <- c("here", "tidyverse", 
                  'emmeans', 'glmmTMB',
                  'patchwork')


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

theme_set(theme_bw())

#get scale_df function
source(here("00_functions",
            "tidy_functions.R"))

# Load data ---------------------------------------------------------------

#need:
#summaries of parameter effects
#raw data
fish_raw <- read.csv(here('01_sbc_fish',
                          'data_outputs',
                          'MSAM',
                          'all_fish_data.csv'))

fish_raw2 <- fish_raw %>%
  distinct(SITE_TRANS, YEAR, MONTH, VIS2, REP, yrID, siteID)

fish_sizes <- read.csv(here('01_sbc_fish',
                            'data_outputs',
                            'MSAM',
                            'all_fish_size_data.csv'))

fish_sum <- readRDS(here('01_sbc_fish',
                         'monsoon',
                         'fish_MSAM',
                         'outputs',
                         'fish_detection_summary.RDS'))


# Get scaled dfs for fish dataset -----------------------------------------

b0 <- as.data.frame(fish_sum$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(parm == "mu.a0") %>%
  dplyr::select(`50%`) %>%
  as_vector()

#### VISIBILITY
vis <- scale_df(x = fish_raw2$VIS2,
                length = 20,
                name = "vis")

bvis <- as.data.frame(fish_sum$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(parm == "a1.Vis") %>%
  dplyr::select(`50%`) %>%
  as_vector()

regvis <- vis %>%
  mutate(reg = b0 + bvis*varS,
         plogis_reg = plogis(reg))

ggplot() +
  geom_line(data = regvis, aes(x = vis, y = plogis_reg)) +
  labs(x = "Dive visibility (m)",
       y = "Detection probability")
  

##### BODYSIZE
fishsz <- scale_df(x = fish_sizes$AVG_SIZE,
                length = 20,
                name = "size")

bfishsz <- as.data.frame(fish_sum$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(parm == "a2.Size") %>%
  dplyr::select(`50%`) %>%
  as_vector()

regfsz <- fishsz %>%
  mutate(reg = b0 + bfishsz*varS,
         plogis_reg = plogis(reg))

ggplot() +
  geom_line(data = regfsz, aes(x = size, y = plogis_reg)) +
  labs(x = "Body size (cm)",
       y = "Detection probability")
