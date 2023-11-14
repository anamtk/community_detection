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
#summaries of parameter effect
fish_sum <- readRDS(here('01_sbc_fish',
                         'monsoon',
                         'fish_MSAM',
                         'outputs',
                         'fish_detection_summary.RDS'))

bird_sum <- readRDS(here('02_konza_birds',
                         'monsoon',
                         'MSAM',
                         'outputs',
                         'bird_detection_summary.RDS'))

plant_sum <- readRDS(here('04_nps_plants',
                          'monsoon',
                          'nps_MSAM',
                          'outputs_yrsite',
                          'nps_detection_summary.RDS'))

# Effect plots ------------------------------------------------------------


# Fish --------------------------------------------------------------------


fish_effects <- as.data.frame(fish_sum$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(parm %in% c("a1.Vis", "a2.Size"))

(fishdetect <- ggplot(fish_effects, aes(x = parm, y = `50%`)) +
  geom_hline(yintercept = 0, linetype = 2, size = 0.75) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), size = 0.75, width = 0) +
  labs(x = "Detection covariate",
       y = "Covariate effect \n (Median and 95% BCI)",
       title = "SBC fish") +
  scale_x_discrete(labels = c("Dive visibility", "Fish size")) + 
  coord_flip() +
  theme(axis.text = element_text(size = 12),
        axis.title= element_text(size = 15),
        plot.title.position = "panel",
        plot.title = element_text(hjust = 0.5)))

# ggsave(plot = fishdetect,
#        filename = here('pictures',
#                        'detection_models',
#                        'detection_covariate_effects.jpg'),
#        height = 4, 
#        width = 6,
#        units = "in")


# Birds -------------------------------------------------------------------



bird_effects <- as.data.frame(bird_sum$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(parm %in% c("a1.Effort", "a2.Size"))

(birddetect <- ggplot(bird_effects, aes(x = parm, y = `50%`)) +
    geom_hline(yintercept = 0, linetype = 2, size = 0.75) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), size = 0.75, width = 0) +
    labs(x = "Detection covariate",
         y = "Covariate effect \n (Median and 95% BCI)",
         title = "KNZ birds") +
    scale_x_discrete(labels = c("Survey effort", "Bird size")) + 
    coord_flip() +
    theme(axis.text = element_text(size = 12),
          axis.title= element_text(size = 15),
          plot.title.position = "panel",
          plot.title = element_text(hjust = 0.5)))


# Plants ------------------------------------------------------------------


plant_effects <- as.data.frame(plant_sum$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(str_detect(parm, "a1.Cover|a2.LifeGroup"))

plant_effects1 <- plant_effects %>%
  filter(parm == "a1.Cover")

(plantdetect1 <- ggplot(plant_effects1, aes(x = parm, y = `50%`)) +
    geom_hline(yintercept = 0, linetype = 2, size = 0.75) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), size = 0.75, width = 0) +
    labs(x = "Detection covariate",
         y = "Covariate effect \n (Median and 95% BCI)",
         title = "PFNP Plants") +
    scale_x_discrete(labels = c("Average species cover")) + 
    coord_flip() +
    theme(axis.text = element_text(size = 12),
          axis.title= element_text(size = 15),
          plot.title.position = "panel",
          plot.title = element_text(hjust = 0.5)))

plant_effects2 <- plant_effects %>%
  filter(str_detect(parm, "a2")) %>%
  filter(parm != "a2.LifeGroup[1]")

(plantdetect2 <- ggplot(plant_effects2, aes(x = parm, y = `50%`)) +
    geom_hline(yintercept = 0, linetype = 2, size = 0.75) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), size = 0.75, width = 0) +
    labs(x = "Life group",
         y = "Covariate effect \n (Median and 95% BCI)") +
    scale_x_discrete(labels = c("Perennial forb", "Perennial grass", 'Perennial shrub',
                                'Annual grass', 'Biennial forb', 'Perennial cacti',
                                'Perennial succulent')) + 
    theme(axis.text = element_text(size = 12),
          axis.title= element_text(size = 15),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title.position = "panel",
          plot.title = element_text(hjust = 0.5)) +
  annotate(geom = "text", x = 3, y = 45, label = "*", size = 8))

plants <- plantdetect1 + plantdetect2

# String them together ----------------------------------------------------

#I NEED HELP MAKING THESE PRETTTTHHHH

fishdetect / birddetect/ plants


# Old code for partial plots ----------------------------------------------



# ggsave(plot = birddetect,
#        filename = here('pictures',
#                        'detection_models',
#                        'detection_covariate_effects.jpg'),
#        height = 4, 
#        width = 6,
#        units = "in")

# Get scaled dfs for fish dataset -----------------------------------------
# 
# #raw data
# fish_raw <- read.csv(here('01_sbc_fish',
#                           'data_outputs',
#                           'MSAM',
#                           'all_fish_data.csv'))
# 
# fish_raw2 <- fish_raw %>%
#   distinct(SITE_TRANS, YEAR, MONTH, VIS2, REP, yrID, siteID)
# 
# fish_sizes <- read.csv(here('01_sbc_fish',
#                             'data_outputs',
#                             'MSAM',
#                             'all_fish_size_data.csv'))
# 
# b0 <- as.data.frame(fish_sum$quantiles) %>%
#   rownames_to_column(var = "parm") %>%
#   filter(parm == "mu.a0") %>%
#   dplyr::select(`50%`) %>%
#   as_vector()
# 
# #### VISIBILITY
# vis <- scale_df(x = fish_raw2$VIS2,
#                 length = 20,
#                 name = "vis")
# 
# bvis <- as.data.frame(fish_sum$quantiles) %>%
#   rownames_to_column(var = "parm") %>%
#   filter(parm == "a1.Vis") %>%
#   dplyr::select(`50%`) %>%
#   as_vector()
# 
# regvis <- vis %>%
#   mutate(reg = b0 + bvis*varS,
#          plogis_reg = plogis(reg))
# 
# (fishvis <- ggplot() +
#     geom_line(data = regvis, aes(x = vis, y = plogis_reg), size = 1) +
#     labs(x = "Dive visibility (m)",
#          y = "Detection probability") + 
#     ylim(0, 0.5) )
# 
# 
# ##### BODYSIZE
# fishsz <- scale_df(x = fish_sizes$AVG_SIZE,
#                    length = 20,
#                    name = "size")
# 
# bfishsz <- as.data.frame(fish_sum$quantiles) %>%
#   rownames_to_column(var = "parm") %>%
#   filter(parm == "a2.Size") %>%
#   dplyr::select(`50%`) %>%
#   as_vector()
# 
# regfsz <- fishsz %>%
#   mutate(reg = b0 + bfishsz*varS,
#          plogis_reg = plogis(reg))
# 
# (fishsize <- ggplot() +
#     geom_line(data = regfsz, aes(x = size, y = plogis_reg), size = 1) +
#     labs(x = "Body size (cm)",
#          y = "Detection probability") +
#     ylim(0, 0.5) +
#     theme(axis.title.y = element_blank()))
# 
# fishvis + fishsize


