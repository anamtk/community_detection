#SAM model results
#Ana Miller-ter Kuile
#October 27, 2023

#this script generates results figures for SAM models

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

theme_set(theme_classic())
theme_update(plot.title.position = "plot",
             panel.grid = element_blank(),
             axis.title = element_text(size = 6),
             axis.text = element_text(size = 5),
             plot.title = element_text(size = 7),
             legend.text = element_text(size = 4),
             legend.key.size = unit(0.25, "cm"),
             legend.key = element_blank(),
             legend.background = element_blank())

source(here('00_functions',
            'tidy_functions.R'))

source(here('00_functions',
            'plot_functions.R'))
# Load data ---------------------------------------------------------------

fish_sam <- readRDS(here('01_sbc_fish',
                         'monsoon',
                         'SAM',
                         'outputs',
                         'fish_SAM_summary.RDS'))

fish_bray <- read.csv(here("01_sbc_fish",
                           "data_outputs",
                           'SAM',
                           'data_prep',
                           "stability_metrics_with_covariates_long.csv"))

sev_sam <- readRDS(here('03_sev_grasshoppers',
                        'monsoon',
                        'SAM',
                        'outputs',
                        'sev_SAM_summary.RDS'))

sev_bray <- read.csv(here("03_sev_grasshoppers",
                          "data_outputs",
                          "SAM",
                          'data_prep',
                          'sev_stability_metrics_with_covariates.csv'))

bird_sam <- readRDS(here('02_konza_birds',
                         'data_outputs',
                         'SAM',
                         'model_outputs',
                         'knz_SAM_summary.RDS'))

bird_bray <- read.csv(here('02_konza_birds',
                           'data_outputs',
                           "SAM",
                           'data_prep',
                           'stability_metrics_with_covariates.csv'))

plant_sam <- readRDS(here('04_nps_plants',
                          'monsoon',
                          'nps_SAM',
                          'outputs',
                          'nps_SAM_summary.RDS'))

plant_diss <- read.csv(here('04_nps_plants',
                            'data_outputs',
                            'SAM',
                            'data_prep',
                            'nps_stability_metrics_with_covariates.csv'))

plant_diss <- plant_diss %>%
  unite(c(Plot, Transect, Quadrat),
        col = "plot_trans_quad",
        sep = "_", 
        remove = F)
# Effect plots ------------------------------------------------------------

fisheffectsplot <- effects_plot_fun(fish_sam) +
  labs(title = "(a) SBC fish") +
  scale_x_continuous(limits = c(-0.5, 0.5)) +
  scale_y_discrete(labels = c("Kelp biomass", "Temperature")) 

birdeffectsplot <- effects_plot_fun(bird_sam) + 
  labs(title = "(b) KNZ birds") +
  scale_x_continuous(limits = c(-1.5, 1.5), breaks = seq(from = -1.5, to = 1.5, by = 0.5)) +
  scale_y_discrete(labels = c("Precipitation", "Temperature"))

seveffectsplot <- effects_plot_fun(sev_sam)  +
  labs(title = "(c) SEV grasshoppers") +
  scale_x_continuous(limits = c(-3.5, 3.5), breaks = seq(from = -3, to = 3, by = 1)) +
  scale_y_discrete(labels = c("Temperature", "Precipitation", "Plant biomass")) 

planteffectsplot <- effects_plot_fun(plant_sam) +
  labs(title = "(d) PFNP plants") +
  scale_x_continuous(limits = c(-2.5, 2.5)) +
  scale_y_discrete(labels = c("Vapor pressure deficit", "Precipitation")) 

effectplots_together <- fisheffectsplot / birdeffectsplot / seveffectsplot / planteffectsplot

# ggsave(plot = effectplots_together,
#        filename = here('pictures',
#                        'sam_models',
#                        'sam_covariate_effects.jpg'),
#        height = 14, 
#        width = 8,
#        units = "cm",
#        dpi = 300)

# Partial plots -----------------------------------------------------------


# Fish partial plots ------------------------------------------------------

fish_title <- ggplot(data.frame(l = "SBC fish", x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), size = 4) + 
  theme_void() +
  coord_cartesian(clip = "off")

#temperature
fisht <- partial_plot_fun(model = fish_sam, 
                          covariate = 'b[2]', 
                          df = fish_bray, 
                          ID = 'SITE_TRANS', 
                          yearID = 'YEAR', 
                          start = 'TEMP_C', 
                          end = 'TEMP_C_l10',
                          weight = "wB",
                          diss = as.name('bray')) +
  labs(x = "Temperature",
       y = "Bray-Curtis Dissimilarity",
       title = "(a)") + 
  theme(plot.title.position = "plot")

###WEIGHTS
fish_tweights <- as.data.frame(fish_sam$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(str_detect(parm, "wB")) %>%
  mutate(season = case_when(parm %in% c("wB[1]", "wB[3]", "wB[5]",
                                        'wB[7]', 'wB[9]','wB[11]') ~ "Warm",
                            parm %in% c("wB[2]", "wB[4]", 'wB[6]',
                                        'wB[8]', 'wB[10]') ~ "Cold")) %>%
  mutate(year = case_when(parm == "wB[1]" ~ 0,
                          parm %in% c("wB[2]", 'wB[3]') ~ 1,
                          parm %in% c('wB[4]', 'wB[5]') ~ 2,
                          parm %in% c("wB[6]", 'wB[7]') ~ 3,
                          parm %in% c("wB[8]", 'wB[9]') ~ 4,
                          parm %in% c("wB[10]", 'wB[11]') ~ 5,
                          TRUE ~ NA_real_)) %>% 
  mutate(above = case_when(
    `50%` > 1/11 ~ "yes",
    TRUE ~ "no"
  )) %>% 
  complete(season, year) %>% 
  unite("color", season, above, sep = "-", remove = FALSE)

warmcol <- '#d8b365'
coldcol <- '#5ab4ac'

warmalpha <- '#f0e1c4'
coldalpha <- '#c2e1df'

fish_tweights_plot <- fish_tweights %>%
  ggplot(aes(x = year, y = `50%`, color = color, shape = color)) +
  geom_hline(yintercept = 1/11, linetype = 2, linewidth = 0.1) +
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`),
                  position = position_dodge(width = 0.5), size = 0.4) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5)) +
  scale_color_manual(values = c("Warm-yes" = warmcol, "Cold-yes" = coldcol,
                                "Warm-no" = warmalpha, "Cold-no" = coldalpha),
                     breaks = c("Warm-yes", "Cold-yes"),
                     labels = c("Warm-yes" = "Warm", "Cold-yes" = "Cold")) +
  scale_shape_manual(values = c("Warm-yes" = 17, "Cold-yes" = 16,
                                "Warm-no" = 17, "Cold-no" = 16),
                     breaks = c("Warm-yes", "Cold-yes"),
                     labels = c("Warm-yes" = "Warm", "Cold-yes" = "Cold")) +
  labs(x = "Years into the past",
       y = "Importance weight",
       shape = "", color = "",
       title = "(b)") +
  theme(legend.position = c(0.92, 1.1))

# Bird partial plots ------------------------------------------------------

bird_title <- ggplot(data.frame(l = "KNZ birds", x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), size = 4) + 
  theme_void() +
  coord_cartesian(clip = "off")

#temperature
birdt <- partial_plot_fun(model = bird_sam, 
                          covariate = 'b[1]', 
                          df = bird_bray, 
                          ID= 'WATERSHED', 
                          yearID = 'RECYEAR', 
                          start = 'TAVE', 
                          end = 'TAVE_l5',
                          weight = "wA",
                          diss = as.name('bray')) +
  labs(x = "Temperature",
       y = "Bray-Curtis Dissimilarity",
       title = "(c)") 

###WEIGHTS
bird_tweights <- as.data.frame(bird_sam$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(str_detect(parm, "wA")) %>%
  mutate(season = case_when(parm %in% c("wA[1]", "wA[3]", "wA[5]") ~ "Cold",
                            parm %in% c("wA[2]", "wA[4]", "wA[6]") ~ "Warm")) %>%
  mutate(year = case_when(parm %in% c("wA[1]", "wA[2]") ~ 0,
                          parm %in% c('wA[3]', 'wA[4]') ~ 1,
                          parm %in% c('wA[5]', 'wA[6]') ~ 2)) %>% 
  mutate(above = case_when(
    `50%` > 1/6 ~ "yes",
    TRUE ~ "no"
  )) %>% 
  complete(season, year) %>% 
  unite("color", season, above, sep = "-", remove = FALSE)

bird_tweights_plot <- bird_tweights %>%
  ggplot(aes(x = year, y = `50%`, color = color, shape = color)) +
  geom_hline(yintercept = 1/6, linetype = 2, linewidth = 0.1) +
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`),
                  position = position_dodge(width = 0.5), size = 0.4) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5)) +
  scale_color_manual(values = c("Warm-yes" = warmcol, "Cold-yes" = coldcol,
                                "Warm-no" = warmalpha, "Cold-no" = coldalpha),
                     breaks = c("Warm-yes", "Cold-yes"),
                     labels = c("Warm-yes" = "Warm", "Cold-yes" = "Cold")) +
  scale_shape_manual(values = c("Warm-yes" = 17, "Cold-yes" = 16,
                                "Warm-no" = 17, "Cold-no" = 16),
                     breaks = c("Warm-yes", "Cold-yes"),
                     labels = c("Warm-yes" = "Warm", "Cold-yes" = "Cold")) +
  labs(x = "Years into the past",
       y = "Importance weight",
       shape = "", color = "",
       title = "(d)") +
  theme(legend.position = c(0.92, 1.1))

#PPT

#ppt
birdp <- partial_plot_fun(model = bird_sam, 
                          covariate = 'b[2]', 
                          df = bird_bray, 
                          ID = 'WATERSHED', 
                          yearID = 'RECYEAR', 
                          start = 'PPT', 
                          end = 'PPT_l5',
                          weight = "wB",
                          diss = as.name('bray')) +
  labs(x = "Precipitation",
       y = "Bray-Curtis Dissimilarity",
       title = "(e)")

###WEIGHTS
bird_pweights <- as.data.frame(bird_sam$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(str_detect(parm, "wB")) %>%
  mutate(season = case_when(parm %in% c("wB[1]", "wB[3]", "wB[5]") ~ "Dry",
                            parm %in% c("wB[2]", "wB[4]", "wB[6]") ~ "Wet")) %>%
  mutate(year = case_when(parm %in% c("wB[1]", "wB[2]") ~ 0,
                          parm %in% c('wB[3]', 'wB[4]') ~ 1,
                          parm %in% c('wB[5]', 'wB[6]') ~ 2)) %>% 
  mutate(above = case_when(
    `50%` > 1/6 ~ "yes",
    TRUE ~ "no"
  )) %>% 
  complete(season, year) %>% 
  unite("color", season, above, sep = "-", remove = FALSE)

drycol <- "#E6922D"
wetcol <- "#465C66"

dryalpha <- "#f5d6b5"
wetalpha <- "#bdc2c5"

bird_pweights_plot <- bird_pweights %>%  
  ggplot(aes(x = year, y = `50%`, color = color, shape = color)) +
  geom_hline(yintercept = 1/6, linetype = 2, linewidth = 0.1) +
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`),
                  position = position_dodge(width = 0.5), size = 0.4) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5)) +
  scale_color_manual(values = c("Dry-yes" = drycol, "Wet-yes" = wetcol,
                                "Dry-no" = dryalpha, "Wet-no" = wetalpha),
                     breaks = c("Dry-yes", "Wet-yes"),
                     labels = c("Dry-yes" = "Dry", "Wet-yes" = "Wet")) +
  scale_shape_manual(values = c("Dry-yes" = 15, "Wet-yes" = 18,
                                "Dry-no" = 15, "Wet-no" = 18),
                     breaks = c("Dry-yes", "Wet-yes"),
                     labels = c("Dry-yes" = "Dry", "Wet-yes" = "Wet")) +
  labs(x = "Years into the past",
       y = "Importance weight",
       shape = "", color = "",
       title = "(f)") +
  theme(legend.position = c(0.92, 1.1))

# Grasshopper partial plots -----------------------------------------------

hopper_title <- ggplot(data.frame(l = "SEV grasshoppers", x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), size = 4) + 
  theme_void() +
  coord_cartesian(clip = "off")

#temperature
hoppert <- partial_plot_fun(model = sev_sam, 
                            covariate = 'b[1]', 
                            df = sev_bray, 
                            ID= 'site_web_trans', 
                            yearID = 'YEAR', 
                            start = 'Temp', 
                            end = 'Temp_l5',
                            weight = "wA",
                            diss = as.name('bray')) +
  labs(x = "Temperature",
       y = "Bray-Curtis Dissimilarity",
       title = "(g)")

###WEIGHTS
hopper_tweights <- as.data.frame(sev_sam$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(str_detect(parm, "wA")) %>%
  mutate(season = case_when(parm %in% c("wA[1]", "wA[3]", "wA[5]") ~ "Cold",
                            parm %in% c("wA[2]", "wA[4]", "wA[6]") ~ "Warm")) %>%
  mutate(year = case_when(parm %in% c("wA[1]", "wA[2]") ~ 0,
                          parm %in% c('wA[3]', 'wA[4]') ~ 1,
                          parm %in% c('wA[5]', 'wA[6]') ~ 2)) %>% 
  mutate(above = case_when(
    `50%` > 1/6 ~ "yes",
    TRUE ~ "no"
  )) %>% 
  complete(season, year) %>% 
  unite("color", season, above, sep = "-", remove = FALSE)

hopper_tweights_plot <- hopper_tweights %>%
  ggplot(aes(x = year, y = `50%`, color = color, shape = color)) +
  geom_hline(yintercept = 1/6, linetype = 2, linewidth = 0.1) +
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`),
                  position = position_dodge(width = 0.5), size = 0.4) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5)) +
  scale_color_manual(values = c("Warm-yes" = warmcol, "Cold-yes" = coldcol,
                                "Warm-no" = warmalpha, "Cold-no" = coldalpha),
                     breaks = c("Warm-yes", "Cold-yes"),
                     labels = c("Warm-yes" = "Warm", "Cold-yes" = "Cold")) +
  scale_shape_manual(values = c("Warm-yes" = 17, "Cold-yes" = 16,
                                "Warm-no" = 17, "Cold-no" = 16),
                     breaks = c("Warm-yes", "Cold-yes"),
                     labels = c("Warm-yes" = "Warm", "Cold-yes" = "Cold")) +
  labs(x = "Years into the past",
       y = "Importance weight",
       shape = "", color = "",
       title = "(h)") +
  theme(legend.position = c(0.92, 1.1))

#PPT

hopperPPT <- partial_plot_fun(model = sev_sam, 
                              covariate = 'b[2]', 
                              df = sev_bray, 
                              ID= 'site_web_trans', 
                              yearID = 'YEAR', 
                              start = 'PPT', 
                              end = 'PPT_l5',
                              weight = "wB",
                              diss = as.name('bray')) +
  labs(x = "Precipitation",
       y = "Bray-Curtis Dissimilarity",
       title = "(i)")

###WEIGHTS
hopper_ppt_weights <- as.data.frame(sev_sam$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(str_detect(parm, "wB")) %>%
  mutate(season = case_when(parm %in% c("wB[1]", "wB[3]", "wB[5]", "wB[7]",
                                        'wB[9]', 'wB[11]') ~ "Wet",
                            parm %in% c("wB[2]", "wB[4]", "wB[6]",
                                        'wB[8]', 'wB[10]') ~ "Dry")) %>%
  mutate(year = case_when(parm %in% c("wB[1]", "wB[2]") ~ 0,
                          parm %in% c('wB[3]', 'wB[4]') ~ 1,
                          parm %in% c('wB[5]', 'wB[6]') ~ 2,
                          parm %in% c('wB[7]', 'wB[8]') ~ 3,
                          parm %in% c('wB[9]', 'wB[10]') ~ 4,
                          parm %in% c('wB[11]') ~ 5)) %>% 
  mutate(above = case_when(
    `50%` > 1/12 ~ "yes",
    TRUE ~ "no"
  )) %>% 
  complete(season, year) %>% 
  unite("color", season, above, sep = "-", remove = FALSE)

hopper_ppt_weights_plot <- hopper_ppt_weights %>%  
  ggplot(aes(x = year, y = `50%`, color = color, shape = color)) +
  geom_hline(yintercept = 1/12, linetype = 2, linewidth = 0.1) +
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`),
                  position = position_dodge(width = 0.5), size = 0.4) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5)) +
  scale_color_manual(values = c("Dry-yes" = drycol, "Wet-yes" = wetcol,
                                "Dry-no" = dryalpha, "Wet-no" = wetalpha),
                     breaks = c("Dry-yes", "Wet-yes"),
                     labels = c("Dry-yes" = "Dry", "Wet-yes" = "Wet")) +
  scale_shape_manual(values = c("Dry-yes" = 15, "Wet-yes" = 18,
                                "Dry-no" = 15, "Wet-no" = 18),
                     breaks = c("Dry-yes", "Wet-yes"),
                     labels = c("Dry-yes" = "Dry", "Wet-yes" = "Wet")) +
  labs(x = "Years into the past",
       y = "Importance weight",
       shape = "", color = "",
       title = "(j)") +
  theme(legend.position = c(0.92, 1.1))

#NPP

# hopperNPP <- partial_plot_fun(model = sev_sam, 
#                               covariate = 'b[3]', 
#                               df = sev_bray, 
#                               ID = 'site_web_trans', 
#                               yearID = 'YEAR', 
#                               start = 'NPP', 
#                               end = 'NPP_l5',
#                               weight = "wC",
#                               diss = as.name('bray')) +
#   labs(x = "Plant biomass",
#        y = "Bray-Curtis Dissimilarity",
#        title = "SEV grasshoppers") + 
#   theme(plot.title.position = "panel",
#         plot.title = element_text(hjust = 0.5))

###WEIGHTS
# hopper_npp_weights <- as.data.frame(sev_sam$quantiles) %>%
#   rownames_to_column(var = "parm") %>%
#   filter(str_detect(parm, "wC")) %>%
#   mutate(season = case_when(parm %in% c("wC[1]", "wC[3]", "wC[5]", "wC[7]",
#                                         'wC[9]', 'wC[11]') ~ "Wet",
#                             parm %in% c("wC[2]", "wC[4]", "wC[6]",
#                                         'wC[8]', 'wC[10]') ~ "Dry")) %>%
#   mutate(year = case_when(parm %in% c("wC[1]", "wC[2]") ~ 0,
#                           parm %in% c('wC[3]', 'wC[4]') ~ 1,
#                           parm %in% c('wC[5]', 'wC[6]') ~ 2,
#                           parm %in% c('wC[7]', 'wC[8]') ~ 3,
#                           parm %in% c('wC[9]', 'wC[10]') ~ 4,
#                           parm %in% c('wC[11]') ~ 5))

# (hopper_npp_weights_plot <- hopper_npp_weights %>%
#     ggplot(aes(x = year, y= `50%`, shape = season)) +
#     geom_hline(yintercept = 1/12, linetype = 2) +
#     geom_point(position = position_dodge(width = 0.5), size = 3) +
#     geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), 
#                   position = position_dodge(width = 0.5),
#                   width = 0) +
#     scale_x_continuous(breaks = c(0, 1,2,3,4,5)) +
#     scale_color_manual(values = c(Warm = warmcol, Cold = coldcol)) +
#     labs(x = "Years into the past",
#          y = "Importance weight\n(median and 95% BCI)"))

# (hoppernppgraphs <- hopperNPP + hopper_npp_weights_plot)

# Plant partial plots -----------------------------------------------------

plant_title <- ggplot(data.frame(l = "PFNP plants", x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), size = 4) + 
  theme_void() +
  coord_cartesian(clip = "off")

#temperature
plantp <- partial_plot_fun(model = plant_sam, 
                           covariate = 'b[2]', 
                           df = plant_diss, 
                           ID= 'plot_trans_quad', 
                           yearID = 'EventYear', 
                           start = 'PPT', 
                           end = 'PPT_l20',
                           weight = "wB",
                           diss = as.name('mean')) +
  labs(x = "Precipitation",
       y = "Community turnover",
       title = "(k)") 

###WEIGHTS

# season == "ppt_monsoon" ~ 1,
# season == "ppt_earlysummer" ~ 2,
# season == "ppt_spring" ~ 3,
# ppt_winter == 4

#21 total weights right now:
plant_pweights <- as.data.frame(plant_sam$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(str_detect(parm, "wB")) %>%
  mutate(season = case_when(parm %in% c("wB[1]", "wB[5]", 'wB[9]', 
                                        'wB[13]', 'wB[17]', 'wB[21]') ~ "Monsoon",
                            parm %in% c("wB[2]", 'wB[6]', 'wB[10]',
                                        'wB[14]', 'wB[18]') ~ "Winter",
                            parm %in% c("wB[3]", 'wB[7]', 'wB[11]',
                                        'wB[15]', 'wB[19]') ~ "Spring",
                            parm %in% c("wB[4]", 'wB[8]', 'wB[12]',
                                        'wB[16]', 'wB[20]') ~ "Early Summer")) %>%
  mutate(year = case_when(parm %in% c("wB[1]", 'wB[2]') ~ 0,
                          parm %in% c('wB[3]', 'wB[4]', "wB[5]", 'wB[6]') ~ 1,
                          parm %in% c( 'wB[7]', 'wB[8]','wB[9]', 'wB[10]') ~ 2,
                          parm %in% c('wB[11]', 'wB[12]', 'wB[13]', 'wB[14]') ~ 3,
                          parm %in% c('wB[15]', 'wB[16]', "wB[17]", 'wB[18]') ~ 4,
                          parm %in% c('wB[19]', 'wB[20]', "wB[21]") ~ 5,
                          TRUE ~ NA_real_)) %>% 
  mutate(above = case_when(
    `50%` > 1/21 ~ "yes",
    TRUE ~ "no"
  )) %>% 
  complete(season, year) %>% 
  unite("color", season, above, sep = "-", remove = FALSE)

escol <- "#014034"
mocol <- "#A68B03"
spcol <- "#F29F05"
wicol <- "#D95204"

esalpha <- "#c6d1cf"
moalpha <- "#e6e0c4"
spalpha <- "#fbe1b4"
wialpha <- "#f5d3c6"

plant_pweights_plot <- plant_pweights %>%
  ggplot(aes(x = year, y = `50%`, color = color, shape = color)) +
  geom_hline(yintercept = 1/21, linetype = 2, linewidth = 0.1) +
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`),
                  position = position_dodge(width = 0.5), size = 0.4) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5)) +
  scale_color_manual(values = c("Early Summer-yes" = escol, 
                                "Monsoon-yes" = mocol,
                                "Spring-yes" = spcol,
                                "Winter-yes" = wicol,
                                "Early Summer-no" = esalpha,
                                "Monsoon-no" = moalpha,
                                "Spring-no" = spalpha,
                                "Winter-no" = wialpha),
                     breaks = c("Early Summer-yes", "Monsoon-yes", 
                                "Spring-yes", "Winter-yes"),
                     labels = c("Early Summer", "Monsoon", 
                                "Spring", "Winter")) +
  scale_shape_manual(values = c("Early Summer-yes" = 15, 
                                "Monsoon-yes" = 16,
                                "Spring-yes" = 17,
                                "Winter-yes" = 18,
                                "Early Summer-no" = 15,
                                "Monsoon-no" = 16,
                                "Spring-no" = 17,
                                "Winter-no" = 18),
                     breaks = c("Early Summer-yes", "Monsoon-yes", 
                                "Spring-yes", "Winter-yes"),
                     labels = c("Early Summer", "Monsoon", 
                                "Spring", "Winter")) +
  labs(x = "Years into the past",
       y = "Importance weight",
       shape = "", color = "",
       title = "(l)") +
  theme(legend.position = c(0.92, 1.1))

# Export ------------------------------------------------------------------

sam_partial_plots <- 
  (fish_title) / (fisht + fish_tweights_plot) / 
  (bird_title) / (birdt + bird_tweights_plot) / (birdp + bird_pweights_plot) /
  (hopper_title) / (hoppert + hopper_tweights_plot) / (hopperPPT + hopper_ppt_weights_plot) /
  (plant_title) / (plantp + plant_pweights_plot) +
  plot_layout(heights = c(1.5, 5, 1.5, 5, 5, 1.5, 5, 5, 1.5, 5))

ggsave(filename = here('pictures',
                       'sam_models',
                       'sam_partial_plots_v2.jpg'),
       sam_partial_plots, 
       height = 24,
       width = 18,
       units = "cm",
       dpi = 300)

# Interaction -------------------------------------------------------------

#this interaction is overfitting the data, so i removed it from the model
# blTK <- as.data.frame(fish_sam$quantiles) %>%
#   rownames_to_column(var = "parm") %>%
#   filter(parm == "b[3]") %>%
#   dplyr::select(`50%`) %>%
#   as_vector() 
# 
# ### temperature
# blT <- as.data.frame(fish_sam$quantiles) %>%
#   rownames_to_column(var = "parm") %>%
#   filter(parm == "b[2]") %>%
#   dplyr::select(`50%`) %>%
#   as_vector()
# 
# b0 <- as.data.frame(fish_sam$quantiles) %>%
#   rownames_to_column(var = "parm") %>%
#   filter(str_detect(parm, "b0")) %>%
#   dplyr::select(`50%`) %>%
#   summarise(b0 = mean(`50%`, na.rm = T)) %>%
#   as_vector()
# 
# #get temparutres on scaled scale
# temp_temp <- fish_bray %>%
#   dplyr::select(SITE_TRANS, YEAR, TEMP_C:TEMP_C_l5) %>% #adjust if needed
#   pivot_longer(TEMP_C:TEMP_C_l5,
#                names_to = "lag",
#                values_to = "temp") %>%
#   mutate(temp = scale(temp)) %>%
#   pivot_wider(names_from = "lag",
#               values_from = "temp") %>%
#   dplyr::select(-SITE_TRANS, -YEAR) %>%
#   as.matrix()
# 
# #make scaled data long format to get mean and sd
# tmaxscale <- fish_bray %>%
#   dplyr::select(SITE_TRANS, YEAR, TEMP_C:TEMP_C_l5) %>% #adjust if needed
#   pivot_longer(TEMP_C:TEMP_C_l5,
#                names_to = "lag",
#                values_to = "temp") 
# 
# #get mean and SD of OG data to back-transform stuff
# mean <- mean(tmaxscale$temp, na.rm = T)
# sd <- sd(tmaxscale$temp, na.rm = T)
# 
# #get weights per month
# t_wt <- as.data.frame(fish_sam$quantiles) %>%
#   rownames_to_column(var = "parameter") %>%
#   filter(str_detect(parameter, "wB")) %>%
#   dplyr::select(`50%`) %>%
#   as_vector()
# 
# 
# #get kelp on scaled scale
# kelp_temp <- fish_bray %>%
#   dplyr::select(SITE_TRANS, YEAR, DRY_GM2:DRY_GM2_l5) %>% #adjust if needed
#   pivot_longer(DRY_GM2:DRY_GM2_l5,
#                names_to = "lag",
#                values_to = "kelp") %>%
#   mutate(kelp = scale(kelp)) %>%
#   pivot_wider(names_from = "lag",
#               values_from = "kelp") %>%
#   dplyr::select(-SITE_TRANS, -YEAR) %>%
#   as.matrix()
# 
# #make scaled data long format to get mean and sd
# kelpscale <- fish_bray %>%
#   dplyr::select(SITE_TRANS, YEAR, DRY_GM2:DRY_GM2_l5) %>% #adjust if needed
#   pivot_longer(DRY_GM2:DRY_GM2_l5,
#                names_to = "lag",
#                values_to = "kelp") 
# 
# #get mean and SD of OG data to back-transform stuff
# meank <- mean(kelpscale$kelp, na.rm = T)
# sdk <- sd(kelpscale$kelp, na.rm = T)
# 
# #get weights per month
# k_wt <- as.data.frame(fish_sam$quantiles) %>%
#   rownames_to_column(var = "parameter") %>%
#   filter(str_detect(parameter, "wA")) %>%
#   dplyr::select(`50%`) %>%
#   as_vector()
# 
# #get tmax dataset
# regT <- fish_bray %>%
#   dplyr::select(SITE_TRANS, YEAR, bray, TEMP_C:TEMP_C_l5)
# 
# #multiply months by their weights
# regT$TAnt <- apply(temp_temp, MARGIN = 1, FUN = function(x){sum(x*t_wt)})
# 
# #revert Tmax to OG data scale
# regT <- regT %>%
#   dplyr::select(TAnt, bray, SITE_TRANS, YEAR) %>%
#   mutate(Temp = TAnt*sd + mean)
# 
# #kelp dataset
# regK <- fish_bray %>%
#   dplyr::select(SITE_TRANS, YEAR, bray, DRY_GM2:DRY_GM2_l5)
# 
# #multiply months by their weights
# regK$KAnt <- apply(kelp_temp, MARGIN = 1, FUN = function(x){sum(x*k_wt)})
# 
# #revert Tmax to OG data scale
# regK <- regK %>%
#   dplyr::select(KAnt, bray, SITE_TRANS, YEAR) %>%
#   mutate(Kelp = KAnt*sd + mean)
# 
# #regression
# regB <- regT %>%
#   left_join(regK, by = c("SITE_TRANS", "YEAR", "bray")) 
# 
# temp2 <- scale_df(x = regB$Temp,
#                   length = 20,
#                   name = "temp")
#   
# kelp2 <- scale_df(x = regB$Kelp,
#                   length = 20,
#                   name = "kelp") %>%
#   rename("varK" = "varS",
#          "levelK" = "level")
# 
# tk <- temp2 %>%
#   cross_join(kelp2) %>%
#   mutate(reg = b0 + blT*varS + blTK*varS*varK,
#          plogisreg = plogis(reg))
# 
# a <- ggplot(tk, aes(x = temp, y = kelp, fill = plogisreg)) +
#   geom_tile() +
#   geom_contour(aes(z = plogisreg), color = "white") +
#   scale_fill_viridis_c() +
#   theme(axis.title = element_blank())
# 
# b <- ggplot(regB, aes(x = Kelp)) +
#   geom_boxplot() +
#   coord_flip() +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank())
# 
# c <- ggplot(regB, aes(x = Temp)) +
#   geom_boxplot() +
#   theme(axis.text.y = element_blank(),
#         axis.ticks.y = element_blank())
# 
# (b + a)/(plot_spacer() + c) +
#   plot_layout(widths = c(1, 3),
#               heights = c(3, 1))



# Summaries for paper -----------------------------------------------------

View(fish_sam$quantiles)
View(bird_sam$quantiles)
View(sev_sam$quantiles)
