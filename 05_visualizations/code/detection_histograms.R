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

bird_detect <- readRDS(here('02_konza_birds',
                            'monsoon',
                            'MSAM',
                            'outputs',
                            'bird_p0_summary.RDS'))

plant_detect <- readRDS(here('04_nps_plants',
                             'monsoon',
                             'nps_MSAM',
                             'outputs_yrsite',
                             'nps_p0_summary.RDS'))

# Prep for plotting -------------------------------------------------------

fish_det2 <- as.data.frame(fish_detect$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(parm != "deviance")%>%
  mutate(dataset = "fish")

hop_det2 <- as.data.frame(hopper_detect$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(parm != "deviance") %>%
  mutate(dataset = "hoppers")

bird_det2 <- as.data.frame(bird_detect$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(parm != "deviance") %>%
  mutate(dataset = "birds")

plant_det2 <- as.data.frame(plant_detect$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(parm != "deviance") %>%
  mutate(dataset = "plants")

detect_df <- fish_det2 %>%
  rbind(hop_det2, bird_det2, plant_det2)

# Create histograms -------------------------------------------------------

labs <- c("SBC fish", "SEV grasshoppers", "KNZ birds", "PFNP Plants")
names(labs) <- c("fish", "hoppers", 'birds', 'plants')

detect_df %>%
  mutate(dataset = factor(dataset, levels = c("hoppers", "fish", 'birds', 'plants'))) %>%
ggplot(aes(x = `50%`)) +
  geom_histogram() +
  facet_grid(dataset~.,
             labeller = labeller(dataset = labs)) +
  labs(x = "Species-level detection probability",
       y = "Count") +
  theme(strip.background = element_rect(fill = "white"))

detect_df %>%
  mutate(dataset = factor(dataset, levels = c("fish", 'birds',"hoppers", 'plants'))) %>%
  ggplot(aes(x = dataset, y = `50%`)) +
  geom_boxplot() +
  labs(x = "Dataset", 
       y = "Species-level detection probability") +
  scale_x_discrete(labels = labs) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(plot = last_plot(),
       filename = here("pictures",
                       "detection_models",
                       "species_detection_probabilities.jpg"),
       height = 4,
       width = 5,
       units = "in")
