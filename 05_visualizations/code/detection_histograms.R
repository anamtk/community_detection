#Distribution of detection probabilities
#Ana Miller-ter Kuile
#October 27, 2023

#this script looks at the distribution of detection probabilities for each
#MSAM/MSOM for each dataset


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


# Load data ---------------------------------------------------------------

#need: 
#summary of species-level detection probabilities

fish_detect <- readRDS(here('01_sbc_fish',
                            'monsoon',
                            'fish_MSAM',
                            'outputs',
                            'fish_p0_summary.RDS'))

hopper_detect <- readRDS(here('03_sev_grasshoppers',
                              'monsoon',
                              'MSAM',
                              'outputs',
                              'grasshopper_p_summary.RDS'))

# Prep for plotting -------------------------------------------------------

fish_det2 <- as.data.frame(fish_detect$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(parm != "deviance")%>%
  mutate(dataset = "fish")

hop_det2 <- as.data.frame(hopper_detect$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(parm != "deviance") %>%
  mutate(dataset = "hoppers")

detect_df <- fish_det2 %>%
  rbind(hop_det2)

# Create histograms -------------------------------------------------------

ggplot(detect_df, aes(x = `50%`)) +
  geom_histogram() +
  facet_grid(dataset~.)

