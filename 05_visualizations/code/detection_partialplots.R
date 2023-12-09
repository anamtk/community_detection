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
  filter(parm %in% c("a1.Vis", "a2.Size")) %>% 
  mutate(parm = case_match(parm,
    "a1.Vis" ~ "Dive visibility",
    "a2.Size" ~ "Fish size"
  ))

fishdetect <- ggplot(fish_effects, aes(x = `50%`, y = parm)) +
  geom_vline(xintercept = 0, linetype = 3, size = 0.75, alpha = 0.4) +
  geom_point(size = 2) +
  geom_errorbar(aes(xmin = `2.5%`, xmax = `97.5%`), size = 0.75, width = 0, alpha = 0.6) +
  scale_x_continuous(limits = c(-3.25, 3.25), breaks = c(-3, -2, -1, 0, 1, 2, 3)) +
  labs(x = "Covariate effect",
       y = "SBC fish",
       title = "(a) SBC fish") +
  # coord_flip() +
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 7),
        axis.title.y = element_blank(),
        # axis.title.y = element_text(angle = 0, vjust = 0.5),
        plot.title.position = "plot",
        plot.title = element_text(size = 8),
        panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank()) +
  geom_hline(yintercept = 0.4)

fishdetect

fish_title <- ggplot(data.frame(l = "SBC fish", x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), size = 4.5) + 
  theme_void() +
  coord_cartesian(clip = "off")

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
  filter(parm %in% c("a1.Effort", "a2.Size")) %>% 
  mutate(parm = case_match(parm,
    "a1.Effort" ~ "Survey effort",
    "a2.Size" ~ "Bird size"
  ))

birddetect <- ggplot(bird_effects, aes(x = `50%`, y = parm)) +
  geom_vline(xintercept = 0, linetype = 3, size = 0.75, alpha = 0.4) +
  geom_point(size = 2) +
  geom_errorbar(aes(xmin = `2.5%`, xmax = `97.5%`), size = 0.75, width = 0, alpha = 0.6) +
  scale_x_continuous(limits = c(-1.6, 1.6), breaks = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5)) +
  labs(x = "Covariate effect",
       y = "",
       title = "(b) KNZ birds") +
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 7),
        axis.title.y = element_blank(),
        # axis.title.y = element_text(angle = 0, vjust = 0.5),
        plot.title.position = "plot",
        plot.title = element_text(size = 8),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank()) +
  geom_hline(yintercept = 0.4)

birddetect

bird_title <- ggplot(data.frame(l = "KNZ birds", x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), size = 4.5) + 
  theme_void() +
  coord_cartesian(clip = "off")


# Plants ------------------------------------------------------------------

plant_effects <- as.data.frame(plant_sum$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(str_detect(parm, "a1.Cover|a2.LifeGroup")) %>% 
  mutate(new_parm = case_match(
    parm,
    "a1.Cover" ~ "Average\nspecies\ncover", 
    "a2.LifeGroup[2]" ~ "Perennial forb",
    "a2.LifeGroup[3]" ~ "Perennial grass",
    "a2.LifeGroup[4]" ~ "Perennial shrub",
    "a2.LifeGroup[5]" ~ "Annual grass",
    "a2.LifeGroup[6]" ~ "Biennial forb",
    "a2.LifeGroup[7]" ~ "Perennial cacti",
    "a2.LifeGroup[8]" ~ "Perennial succulent"
  ),
  new_parm = fct_inorder(new_parm)) 

plant_effects1 <- plant_effects %>%
  filter(parm == "a1.Cover")

plantdetect1 <- ggplot(plant_effects1, aes(x = `50%`, y = new_parm)) +
  geom_vline(xintercept = 0, linetype = 3, size = 0.75, alpha = 0.4) +
  geom_point(size = 2) +
  geom_errorbar(aes(xmin = `2.5%`, xmax = `97.5%`), size = 0.75, width = 0) +
  scale_x_continuous(limits = c(-1.75, 1.75), breaks = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5)) +
  labs(x = "Covariate effect",
       y = "",
       title = "(c) PFNP plants: average species cover") +
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 7),
        axis.title.y = element_blank(),
        # axis.title.y = element_text(angle = 0, vjust = 0.5),
        plot.title.position = "plot",
        plot.title = element_text(size = 8),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank()) +
  geom_hline(yintercept = 0.4)

plantdetect1

plant_effects2 <- plant_effects %>%
  filter(str_detect(parm, "a2")) %>%
  filter(parm != "a2.LifeGroup[1]")

plantdetect2 <- ggplot(plant_effects2, aes(x = `50%`, y = new_parm)) +
  geom_vline(xintercept = 0, linetype = 3, size = 0.75, alpha = 0.4) +
  geom_point(size = 2) +
  geom_errorbar(aes(xmin = `2.5%`, xmax = `97.5%`), size = 0.75, width = 0) +
  scale_x_continuous(limits = c(-10, 40)) +
  labs(x = "Covariate effect",
       y = "",
       title = "(d) PFNP plants: Life group") +
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 7),
        axis.title.y = element_blank(),
        # axis.title.y = element_text(angle = 0, vjust = 0.5),
        plot.title.position = "plot",
        plot.title = element_text(size = 8),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank()) +
  geom_hline(yintercept = 0.4)

plantdetect2

plant_title <- ggplot(data.frame(l = "PFNP plants", x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), size = 6) + 
  theme_void() +
  coord_cartesian(clip = "off")

# String them together ----------------------------------------------------

detect_together <- (fishdetect)/(birddetect)/(plantdetect1)/(plantdetect2) +
  plot_layout(heights = c(10, 10, 10, 10))

detect_together

ggsave(plot = detect_together,
       filename = here("pictures",
                       "detection_models",
                       "detection_covariate_effects.jpg"),
       height = 14,
       width = 10,
       units = "cm",
       dpi = 300)

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


