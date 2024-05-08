#SAM model results version 2
#Ana Miller-ter Kuile
#April 2024

#this script generates results figures for SAM models
#with both modeled and raw datasets

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

#fish datasets:
fish_sam <- readRDS(here('01_sbc_fish',
                         'monsoon',
                         'SAM',
                         'modeled',
                         'outputs',
                         'fish_SAM_summary.RDS'))

fish_sam_raw <- readRDS(here("01_sbc_fish",
                             "data_outputs",
                             "SAM",
                             "model_outputs",
                             "fish_SAM_summary_raw.RDS"))

fish_bray <- read.csv(here("01_sbc_fish",
                           "data_outputs",
                           'SAM',
                           'data_prep',
                           "stability_metrics_with_covariates_long.csv"))

#bird datasets:
bird_sam <- readRDS(here('02_konza_birds',
                         'data_outputs',
                         'SAM',
                         'model_outputs',
                         'knz_SAM_summary.RDS'))

bird_sam_raw <- readRDS(here('02_konza_birds',
                             'data_outputs',
                             'SAM',
                             'model_outputs',
                             'knz_SAM_summary_raw.RDS'))

bird_bray <- read.csv(here('02_konza_birds',
                           'data_outputs',
                           "SAM",
                           'data_prep',
                           'stability_metrics_with_covariates.csv'))

#grasshopper datasets:
sev_sam <- readRDS(here('03_sev_grasshoppers',
                        'monsoon',
                        'SAM',
                        'outputs',
                        'sev_SAM_summary.RDS'))

sev_sam_raw <- readRDS(here('03_sev_grasshoppers',
                            'monsoon',
                            'SAM',
                            'raw',
                            'outputs',
                            'sev_SAM_raw_summary.RDS'))

sev_bray <- read.csv(here("03_sev_grasshoppers",
                          "data_outputs",
                          "SAM",
                          'data_prep',
                          'sev_stability_metrics_with_covariates.csv'))

#plant datasets:
plant_sam <- readRDS(here('04_nps_plants',
                          'monsoon',
                          'nps_SAM',
                          'outputs',
                          'nps_SAM_summary.RDS'))

plant_sam_raw <- readRDS(here('04_nps_plants',
                              'data_outputs',
                              'SAM',
                              'model_outputs',
                              'nps_SAM_summary_raw.RDS'))

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

#option for a main figure of effects plots per model

sam_plot_fun <- function(model_model, model_raw){
  
  beta_m <- as.data.frame(model_model$quantiles) %>%
    rownames_to_column(var = "parm") %>%
    filter(str_detect(parm, "b")) %>%
    filter(!str_detect(parm, "b0")) %>%
    filter(!str_detect(parm, 'sig.web')) %>%
    mutate(type = "modeled")
  
  beta_r <- as.data.frame(model_raw$quantiles) %>%
    rownames_to_column(var = "parm") %>%
    filter(str_detect(parm, "b")) %>%
    filter(!str_detect(parm, "b0")) %>%
    filter(!str_detect(parm, 'sig.web')) %>%
    mutate(type = "observed")
  
  all_beta <- beta_m %>% rbind(beta_r)
  
  plot <- ggplot(all_beta, aes(x = `50%`, y = parm, color = type)) +
    geom_vline(xintercept = 0, linetype = 2, alpha = 0.4) +
    geom_point(position=position_dodge(width=0.5)) +
    geom_errorbar(aes(xmin = `2.5%`, xmax = `97.5%`), width = 0,
                  position=position_dodge(width=0.5)) +
    scale_color_manual(values = c(modeled_col, observed_col)) +
    labs(x = "Covariate effect", y = "") 
  
  return(plot)
  
}

# fish --------------------------------------------------------------------

fish_plot <- sam_plot_fun(model_model = fish_sam, model_raw = fish_sam_raw)  +
  scale_y_discrete(labels = c("Kelp biomass", "Temperature")) + 
  labs(title = "SBC fish")

# birds -------------------------------------------------------------------

bird_plot <- sam_plot_fun(model_model = bird_sam, model_raw = bird_sam_raw) +
  scale_y_discrete(labels = c("Precipitation", "Temperature")) + 
  labs(title = "KNZ birds")

# grasshoppers ------------------------------------------------------------

hop_plot <- sam_plot_fun(model_model = sev_sam, model_raw = sev_sam_raw) +
  scale_y_discrete(labels = c("Temperature","Precipitation", "NPP")) + 
  labs(title = "SEV grasshoppers")

# Plant -------------------------------------------------------------------

plant_plot <- sam_plot_fun(model_model = plant_sam, model_raw = plant_sam_raw) +
  scale_y_discrete(labels = c("VPD","Precipitation")) + 
  labs(title = "PFNP plants")


# combine -----------------------------------------------------------------

(fish_plot + bird_plot) / (hop_plot + plant_plot) +
  plot_layout(guides = "collect")

# Weights plots -----------------------------------------------------------

fish_tweights_raw <- as.data.frame(fish_sam_raw$quantiles) %>%
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
  mutate(season_num = case_when(parm == "wB[1]" ~ 1,
                                parm == "wB[2]" ~ 2, 
                                parm == 'wB[3]' ~ 3,
                                parm == 'wB[4]' ~ 4,
                                parm == 'wB[5]' ~ 5,
                                parm == "wB[6]" ~ 6,
                                parm == 'wB[7]' ~ 7,
                                parm == "wB[8]" ~ 8,
                                parm == 'wB[9]' ~ 9,
                                parm == "wB[10]" ~ 10,
                                parm == 'wB[11]' ~ 11,
                                TRUE ~ NA_real_)) %>% 
  mutate(above = case_when(
    `50%` > 1/11 ~ "yes",
    TRUE ~ "no"
  )) %>% 
  complete(season, year) %>% 
  unite("color", season, above, sep = "-", remove = FALSE) %>%
  mutate(type = "raw")

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
  mutate(season_num = case_when(parm == "wB[1]" ~ 1,
                                parm == "wB[2]" ~ 2, 
                                parm == 'wB[3]' ~ 3,
                                parm == 'wB[4]' ~ 4,
                                parm == 'wB[5]' ~ 5,
                                parm == "wB[6]" ~ 6,
                                parm == 'wB[7]' ~ 7,
                                parm == "wB[8]" ~ 8,
                                parm == 'wB[9]' ~ 9,
                                parm == "wB[10]" ~ 10,
                                parm == 'wB[11]' ~ 11,
                          TRUE ~ NA_real_)) %>% 
  mutate(above = case_when(
    `50%` > 1/11 ~ "yes",
    TRUE ~ "no"
  )) %>% 
  complete(season, year) %>% 
  unite("color", season, above, sep = "-", remove = FALSE)%>%
  mutate(type = "modeled")

all_fish_twt <- fish_tweights %>%
  rbind(fish_tweights_raw)

ggplot(all_fish_twt) +
  geom_rect(aes(xmin = 1.5, xmax = 2.5, ymin = 0, ymax = 0.55), 
            fill = "gray92") +
  geom_rect(aes(xmin = 3.5, xmax = 4.5, ymin = 0, ymax = 0.55), 
            fill = "gray92") +
  geom_rect(aes(xmin = 5.5, xmax = 6.5, ymin = 0, ymax = 0.55), 
            fill = "gray92") +
  geom_rect(aes(xmin = 7.5, xmax = 8.5, ymin = 0, ymax = 0.55), 
            fill = "gray92") +
  geom_rect(aes(xmin = 9.5, xmax = 10.5, ymin = 0, ymax = 0.55), 
            fill = "gray92") +
  geom_hline(yintercept = 1/11, linetype = 2, linewidth = 0.1) +
  geom_pointrange(aes(x = season_num, y = `50%`, 
                      color = type, 
                      ymin = `2.5%`, ymax = `97.5%`),
                  position = position_dodge(width = 1), size = 0.4) +
  scale_x_continuous(breaks = c(1, 2, 4,  6,  8,  10), 
                     labels = c(0, 1, 2,  3,  4, 5 )) +
  scale_color_manual(values = c(modeled_col, observed_col)) +
  labs(x = "Years into past",
      y = "Importance weight \n (median and 95% BCI)")



# Effect plots ------------------------------------------------------------

# betas <- as.data.frame(fish_sam$quantiles) %>%
#   rownames_to_column(var = "parm") %>%
#   filter(str_detect(parm, "b")) %>%
#   filter(!str_detect(parm, "b0")) %>%
#   filter(!str_detect(parm, 'sig.web')) %>%
#   mutate(type = "modeled")
# 
# betas2 <- as.data.frame(fish_sam_raw$quantiles) %>%
#   rownames_to_column(var = "parm") %>%
#   filter(str_detect(parm, "b")) %>%
#   filter(!str_detect(parm, "b0")) %>%
#   filter(!str_detect(parm, 'sig.web')) %>%
#   mutate(type = "observed")
# 
# all_betas <- betas %>% rbind(betas2)
# 
# ggplot(all_betas, aes(x = `50%`, y = parm, color = type)) +
#   geom_vline(xintercept = 0, linetype = 2, alpha = 0.4) +
#   geom_point(position=position_dodge(width=0.5)) +
#   geom_errorbar(aes(xmin = `2.5%`, xmax = `97.5%`), width = 0,
#                 position=position_dodge(width=0.5)) +
#   scale_color_manual(values = c(modeled_col, observed_col)) +
#   labs(x = "Covariate effect", y = "") +
#   scale_x_continuous(limits = c(-0.5, 0.5)) +
#   scale_y_discrete(labels = c("Kelp biomass", "Temperature")) 
# 
# 
# 
# # birds -------------------------------------------------------------------
# 
# betas3 <- as.data.frame(bird_sam$quantiles) %>%
#   rownames_to_column(var = "parm") %>%
#   filter(str_detect(parm, "b")) %>%
#   filter(!str_detect(parm, "b0")) %>%
#   filter(!str_detect(parm, 'sig.web')) %>%
#   mutate(type = "modeled")
# 
# betas4 <- as.data.frame(bird_sam_raw$quantiles) %>%
#   rownames_to_column(var = "parm") %>%
#   filter(str_detect(parm, "b")) %>%
#   filter(!str_detect(parm, "b0")) %>%
#   filter(!str_detect(parm, 'sig.web')) %>%
#   mutate(type = "observed")
# 
# all_betas2 <- betas3 %>% rbind(betas4)
# 
# ggplot(all_betas2, aes(x = `50%`, y = parm, color = type)) +
#   geom_vline(xintercept = 0, linetype = 2, alpha = 0.4) +
#   geom_point(position=position_dodge(width=0.5)) +
#   geom_errorbar(aes(xmin = `2.5%`, xmax = `97.5%`), width = 0,
#                 position=position_dodge(width=0.5)) +
#   scale_color_manual(values = c(modeled_col, observed_col)) +
#   labs(x = "Covariate effect", y = "") +
#   #scale_x_continuous(limits = c(-0.5, 0.5)) +
#   scale_y_discrete(labels = c("Precipitation", "Temperature")) 
# 
# # make dataframes ---------------------------------------------------------
# 
# 
# a <- partial_df_fun(model = fish_sam, 
#                     covariate = 'b[2]', 
#                     df = fish_bray, 
#                     ID = 'SITE_TRANS', 
#                     yearID = 'YEAR', 
#                     start = 'TEMP_C', 
#                     end = 'TEMP_C_l10',
#                     weight = "wB",
#                     diss = as.name('bray'))
# 
# b <- partial_df_fun(model = fish_sam_raw, 
#                     covariate = 'b[2]', 
#                     df = fish_bray, 
#                     ID = 'SITE_TRANS', 
#                     yearID = 'YEAR', 
#                     start = 'TEMP_C', 
#                     end = 'TEMP_C_l10',
#                     weight = "wB",
#                     diss = as.name('observed_all'))
# 
# modeled_col <- "#E88C23"
# observed_col <- "#438AA8"
# 
# #Plot regression lines of parital plots
# ggplot() +
#   geom_point(data = a, aes(x = Var, y = bray), 
#              alpha = 0.2, shape = 1,
#              position = position_jitter(),
#              color = modeled_col) +
#   geom_point(data = b, aes(x = Var, y = observed_all),
#              alpha = 0.2, shape = 1,
#              position = position_jitter(),
#              color = observed_col) +
#   geom_line(data = a, aes(x = Var, y = plogisreg), linewidth = 1,
#             color = modeled_col) +
#   geom_line(data = b, aes(x = Var, y = plogisreg), linewidth = 1,
#             color = observed_col)+
#   theme(panel.grid = element_blank(),
#         plot.title.position = "plot") +
#   labs(x = "Temperature",
#        y = "Bray-Curtis Dissimilarity",
#        title = "(a)") 
# 
# c <- partial_df_fun(model = fish_sam, 
#                     covariate = 'b[1]', 
#                     df = fish_bray, 
#                     ID = 'SITE_TRANS', 
#                     yearID = 'YEAR', 
#                     start = 'DRY_GM2', 
#                     end = 'DRY_GM2_l5',
#                     weight = "wA",
#                     diss = as.name('bray'))
# 
# d <- partial_df_fun(model = fish_sam_raw, 
#                     covariate = 'b[1]', 
#                     df = fish_bray, 
#                     ID = 'SITE_TRANS', 
#                     yearID = 'YEAR', 
#                     start = 'DRY_GM2', 
#                     end = 'DRY_GM2_l5',
#                     weight = "wA",
#                     diss = as.name('observed_all'))
# 
# ggplot() +
#   geom_point(data = c, aes(x = Var, y = bray), 
#              alpha = 0.2, shape = 1,
#              position = position_jitter(),
#              color = modeled_col) +
#   geom_point(data = d, aes(x = Var, y = observed_all),
#              alpha = 0.2, shape = 1,
#              position = position_jitter(),
#              color = observed_col) +
#   geom_line(data = c, aes(x = Var, y = plogisreg), linewidth = 1,
#             color = modeled_col) +
#   geom_line(data = d, aes(x = Var, y = plogisreg), linewidth = 1,
#             color = observed_col)+
#   theme(panel.grid = element_blank(),
#         plot.title.position = "plot") +
#   labs(x = "Kelp Biomass",
#        y = "Bray-Curtis Dissimilarity",
#        title = "(a)") 
# 
