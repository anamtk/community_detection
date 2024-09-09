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
             axis.title = element_text(size = 11),
             axis.text = element_text(size = 10),
             plot.title = element_text(size = 12),
             legend.text = element_text(size = 9),
             legend.key.size = unit(0.25, "cm"),
             legend.key = element_blank(), 
             legend.title=element_blank(),
             legend.background = element_blank())

source(here('00_functions',
            'tidy_functions.R'))

source(here('00_functions',
            'plot_functions.R'))

modeled_col <- "#E88C23"
observed_col <- "#438AA8"
# Load data ---------------------------------------------------------------

#fish datasets:
fish_sam <- readRDS(here("01_sbc_fish",
                         "data_outputs",
                         "SAM",
                         "model_outputs",
                         "fish_SAM_summary.RDS"))

fish_sam_raw <- readRDS(here("01_sbc_fish",
                             "data_outputs",
                             "SAM",
                             "model_outputs",
                             "fish_SAM_summary_raw.RDS"))

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

#cummulative temp and ppt from both models
sev_sam_temp <- readRDS(here('03_sev_grasshoppers',
                        'monsoon',
                        'SAM',
                        'outputs',
                        'sev_SAM_cummtemp.RDS'))

sev_sam_ppt <- readRDS(here('03_sev_grasshoppers',
                             'monsoon',
                             'SAM',
                             'outputs',
                             'sev_SAM_cummppt.RDS'))

sev_sam_npp <- readRDS(here('03_sev_grasshoppers',
                            'monsoon',
                            'SAM',
                            'outputs',
                            'sev_SAM_cummnpp.RDS'))

sev_sam_temp_raw <- readRDS(here('03_sev_grasshoppers',
                             'monsoon',
                             'SAM',
                             'raw',
                             'outputs',
                             'sev_SAM_cummtemp_raw.RDS'))

sev_sam_ppt_raw <- readRDS(here('03_sev_grasshoppers',
                             'monsoon',
                             'SAM',
                             'raw',
                             'outputs',
                             'sev_SAM_cummppt_raw.RDS'))

sev_sam_npp_raw <- readRDS(here('03_sev_grasshoppers',
                            'monsoon',
                            'SAM',
                            'raw',
                            'outputs',
                            'sev_SAM_cummnpp_raw.RDS'))


#plant datasets:
plant_sam <- readRDS(here('04_nps_plants',
                          'data_outputs',
                          'SAM',
                          'model_outputs',
                          'nps_SAM_summary.RDS'))

plant_sam_raw <- readRDS(here('04_nps_plants',
                              'data_outputs',
                              'SAM',
                              'model_outputs',
                              'nps_SAM_summary_raw.RDS'))


# Main intercept plot -----------------------------------------------------

#look at intercept in these datasets as the 
#"expected change" under average covariates
#are the intercepts really different? 

b0_fun <- function(model){
  betas <- as.data.frame(model$quantiles) %>%
    rownames_to_column(var = "parm") %>%
    filter(parm == "b0") %>%
    rename("LCI" = `2.5%`,
           "median" = `50%`,
           'UCI' = `97.5%`) %>%
    dplyr::select(LCI, median, UCI) %>%
    mutate(LCI = plogis(LCI),
           median = plogis(median),
           UCI = plogis(UCI))
  
  return(betas)
}

model_list <- list(fish_sam, fish_sam_raw,
                   bird_sam, bird_sam_raw,
                   sev_sam, sev_sam_raw,
                   plant_sam, plant_sam_raw)

names <- c("SBC fish_modeled", "SBC fish_observed",
           "KNZ birds_modeled", "KNZ birds_observed",
           'SEV grasshoppers_modeled', 
           "SEV grasshoppers_observed",
           'PFNP plants_modeled', "PFNP plants_observed")

names(model_list) <- names

b0_list <- lapply(model_list, b0_fun)

b0_df <- bind_rows(b0_list, .id = "id") %>%
  separate(id, 
           into = c('dataset', 'type'),
           sep = "_") %>%
  mutate(dataset = factor(dataset,
                          levels = c("SBC fish", 
                                     "KNZ birds", 
                                     "SEV grasshoppers", 
                                     "PFNP plants"))) %>%
  mutate(datasetID = case_when(dataset == "SBC fish" ~ "Fish",
                               dataset == "KNZ birds" ~ "Birds",
                               dataset == "SEV grasshoppers" ~ "Grasshoppers",
                               dataset == "PFNP plants" ~ "Plants"))

(intercept_plot <- ggplot(b0_df) +
  geom_pointrange(aes(x = datasetID,
                      y = 1-median,
                      ymin = 1-LCI,
                      ymax = 1-UCI,
                      color = type),
                  position=position_dodge(width=0.5),
                  size = 0.4) +
  scale_color_manual(values = c(modeled_col, observed_col),
                     labels = c("Modeled", "Empirical \n (maximum)")) +
  labs(x = "Dataset", y = "Baseline stability (1-intercept) \n (median and 95% BCI)"))
# Effect plots ------------------------------------------------------------

#option for a main figure of effects plots per model

#TO-DO:
#make a boxplot of uncertainty (modeled
#communities have lower uncertainties)

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
  
  all_beta <- beta_m %>% rbind(beta_r) %>%
    mutate(sig = case_when(`2.5%` < 0 & `97.5%` > 0 ~ 'nonsig',
                            TRUE ~ "sig")) %>%
    mutate(alpha = case_when(sig == "sig" ~ 1,
                             sig == "nonsig" ~ 0.2))
  
  plot <- ggplot(all_beta, aes(x = `50%`, y = parm, color = type)) +
    geom_vline(xintercept = 0, linetype = 2, alpha = 0.4) +
    geom_pointrange(aes(x = `50%`,
                        y = parm, 
                        xmin = `2.5%`,
                        xmax = `97.5%`,
                        color = type, 
                        shape = sig), 
                    show.legend = F,
                    position=position_dodge(width=0.5),
                    size = 0.4) +
    scale_color_manual(values = c(modeled_col, observed_col),
                       labels = c("Modeled", "Empirical \n (maximum)")) +
    scale_shape_manual(values = c('nonsig' = 1, 'sig' = 19)) +
    labs(x = "Covariate effect \n (Median and 95% BCI)", y = "") 
  
  return(plot)
  
}

# fish --------------------------------------------------------------------

(fish_plot <- sam_plot_fun(model_model = fish_sam, model_raw = fish_sam_raw)  +
  scale_y_discrete(labels = c("Kelp biomass", "Temperature")) + 
  labs(title = "Fish"))

# birds -------------------------------------------------------------------

bird_plot <- sam_plot_fun(model_model = bird_sam, model_raw = bird_sam_raw) +
  scale_y_discrete(labels = c("Temperature", "Precipitation")) + 
  labs(title = "Birds")

# grasshoppers ------------------------------------------------------------

hop_plot <- sam_plot_fun(model_model = sev_sam, model_raw = sev_sam_raw) +
  scale_y_discrete(labels = c("Temperature","Precipitation", "Plant biomass")) + 
  labs(title = "Grasshoppers")

# Plant -------------------------------------------------------------------

plant_plot <- sam_plot_fun(model_model = plant_sam, model_raw = plant_sam_raw) +
  scale_y_discrete(labels = c("VPD","Precipitation")) + 
  labs(title = "Plants") 

# combine -----------------------------------------------------------------

(effects_plots_all <- intercept_plot / ((fish_plot + bird_plot) / (hop_plot + plant_plot)) + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'a',
                  tag_prefix = "(",
                  tag_suffix = ")"))

ggsave(filename = here('pictures',
                       'sam_models',
                       'sam_effects_plots.jpg'),
       effects_plots_all, 
       height = 20,
       width = 20,
       units = "cm",
       dpi = 300)
# Weights plots -----------------------------------------------------------

# weight function ---------------------------------------------------------

#Do cumulative weights like figure in OG legacy paper
# grey region = raw data
# points are the modeled data results

#need to do this in a different way
#extract 1000 coda samples of the weights and then
#get a cumulative weight for each iteration and THEN
#do the cumulative weights median and credible interval calculation

#modeled data model
#raw data model
#cummulative weight variable ID
#dataset ID
#variable name
weights_fun <- function(model_m, model_r, weightID, 
                        varID, datasetID){
  
  df <- as.data.frame(model_m$quantiles) %>%
    rownames_to_column(var = "parm") %>%
    filter(str_detect(parm, weightID)) %>%
    separate(parm, into = c("parm", 'lag'),
             sep = "\\[") %>%
    mutate(lag = substr(lag, 1, (nchar(lag)-1))) %>%
    mutate(lag = as.numeric(lag)) %>%
    mutate(variable = varID) %>%
    mutate(dataset = datasetID) %>%
    rename(LCI = `2.5%`,
           median = `50%`,
           UCI = `97.5%`) %>%
    dplyr::select(lag, variable, dataset,
                  LCI, median, UCI)
  
  df2 <- as.data.frame(model_r$quantiles) %>%
    rownames_to_column(var = "parm") %>%
    filter(str_detect(parm, weightID)) %>%
    separate(parm, into = c("parm", 'lag'),
             sep = "\\[") %>%
    mutate(lag = substr(lag, 1, (nchar(lag)-1))) %>%
    mutate(lag = as.numeric(lag)) %>%
    mutate(variable = varID) %>%
    mutate(dataset = datasetID) %>%
    rename(LCI_obs = `2.5%`,
           median_obs = `50%`,
           UCI_obs = `97.5%`) %>%
    dplyr::select(lag, variable, dataset,
                  LCI_obs, median_obs, UCI_obs)
  
  df_all <- df %>%
    left_join(df2, by = c("lag", "variable", 'dataset'))
  
  plot <- ggplot(df_all) +
    geom_ribbon(aes(x = lag, ymin = LCI_obs,
                    ymax = UCI_obs),
                fill = observed_col, alpha = 0.5) +
    geom_line(aes(x = lag, y = median_obs), 
              color = observed_col, alpha = 0.5) +
    geom_pointrange(aes(x = lag, y = median,
                        ymin = LCI, 
                        ymax = UCI),
                    color = modeled_col) +    
    geom_line(aes(x = lag, y = median), 
              color = modeled_col, alpha = 0.5,
              linewidth = 0.6) +
    labs(x = "Seasons into the past", 
         y = "Cumulative seasonal weights \n posterior median and 95% BCI")
  
  return(plot)

}

(fishwts_temp <- weights_fun(model_m = fish_sam,
                             model_r = fish_sam_raw,
                             weightID = "cumm.tempwt",
                             #nlag = 6,
                             var = "Temperature",
                             datasetID = "SBC fish") +
    labs(title = "Fish:\ntemperature") +
    scale_x_continuous(limits = c(1,11), 
                       breaks = seq(1, 11, by = 1),
                       labels = c(0:10)))

sevwts_ppt <- sev_sam_ppt %>%
  left_join(sev_sam_ppt_raw, by = c("lag", "variable",
                                    'dataset'))

(sevwts_pptp <- ggplot(sevwts_ppt) +
  geom_ribbon(aes(x = lag, ymin = LCI_obs,
                  ymax = UCI_obs),
              fill = observed_col, alpha = 0.5) +
  geom_line(aes(x = lag, y = median_obs), 
            color = observed_col, alpha = 0.5) +
  geom_pointrange(aes(x = lag, y = median,
                      ymin = LCI, 
                      ymax = UCI),
                  color = modeled_col) +    
  geom_line(aes(x = lag, y = median), 
            color = modeled_col, alpha = 0.5,
            linewidth = 0.6) +
  ylim(0, 1.01) +
  labs(x = "Seasons into the past", 
       y = "Cumulative seasonal weights \n posterior median and 95% BCI") +
    labs(title = "Grasshoppers:\nprecipitation") +
    scale_x_continuous(limits = c(1,6), 
                       breaks = seq(1, 6, by = 1),
                       labels = c(0:5)) +
    theme(axis.title.y = element_blank()))

sevwts_temp <- sev_sam_temp %>%
  left_join(sev_sam_temp_raw, by = c("lag", "variable",
                                    'dataset'))

(sevwts_tempp <- ggplot(sevwts_temp) +
    geom_ribbon(aes(x = lag, ymin = LCI_obs,
                    ymax = UCI_obs),
                fill = observed_col, alpha = 0.5) +
    geom_line(aes(x = lag, y = median_obs), 
              color = observed_col, alpha = 0.5) +
    geom_pointrange(aes(x = lag, y = median,
                        ymin = LCI, 
                        ymax = UCI),
                    color = modeled_col) +    
    geom_line(aes(x = lag, y = median), 
              color = modeled_col, alpha = 0.5,
              linewidth = 0.6) +
    labs(x = "Seasons into the past", 
         y = "Cumulative seasonal weights \n posterior median and 95% BCI") +
    labs(title = "Grasshoppers:\ntemperature") +
    scale_x_continuous(limits = c(1,6), 
                       breaks = seq(1, 6, by = 1),
                       labels = c(0:5))) +
  theme(axis.title.y = element_blank())

sevwts_npp <- sev_sam_npp %>%
  left_join(sev_sam_npp_raw, by = c("lag", "variable",
                                     'dataset'))

(sevwts_nppp <- ggplot(sevwts_npp) +
    geom_ribbon(aes(x = lag, ymin = LCI_obs,
                    ymax = UCI_obs),
                fill = observed_col, alpha = 0.5) +
    geom_line(aes(x = lag, y = median_obs), 
              color = observed_col, alpha = 0.5) +
    geom_pointrange(aes(x = lag, y = median,
                        ymin = LCI, 
                        ymax = UCI),
                    color = modeled_col) +    
    geom_line(aes(x = lag, y = median), 
              color = modeled_col, alpha = 0.5,
              linewidth = 0.6) +
    labs(x = "Seasons into the past", 
         y = "Cumulative seasonal weights \n posterior median and 95% BCI") +
    labs(title = "Grasshoppers:\nPlant biomass") +
    scale_x_continuous(limits = c(1,11), 
                       breaks = seq(1, 11, by = 1),
                       labels = c(0:10))) +
  theme(axis.title.y = element_blank())


fishwts_temp +  sevwts_pptp + sevwts_tempp +
  sevwts_nppp +
  plot_annotation(tag_levels = 'a',
                  tag_prefix = "(",
                  tag_suffix = ")")

ggsave(filename = here('pictures',
                       'sam_models',
                       'sam_weight_plots.jpg'),
       height = 20,
       width = 20,
       units = "cm",
       dpi = 300)
# Old ---------------------------------------------------------------------


# weights_fun3 <- function(model_m, model_r,
#                          weightID, nlag, var, datasetID){
#   
#   df <- as.data.frame(model_m$quantiles) %>%
#     rownames_to_column(var = "parm") %>%
#     filter(str_detect(parm, weightID)) %>%
#     mutate(lag = substr(parm, 4, (nchar(parm) - 1))) %>%
#     mutate(lag = as.numeric(lag)) %>%
#     mutate(variable = var) %>%
#     mutate(dataset = datasetID) %>%
#     rename(LCI = `2.5%`,
#            median = `50%`,
#            UCI = `97.5%`) %>%
#     dplyr::select(lag, variable, dataset,
#                   LCI, median, UCI) %>% 
#     mutate(median_c = accumulate(median, `+`)) %>%
#     rowwise() %>%
#     mutate(LCI_diff = median_c -(median - LCI),
#            UCI_diff = median_c + (UCI - median)) %>%
#     mutate(UCI_diff = case_when(UCI_diff > 1 ~ 1,
#                                 TRUE ~ UCI_diff))
#   
#   
#   df2 <- as.data.frame(model_r$quantiles) %>%
#     rownames_to_column(var = "parm") %>%
#     filter(str_detect(parm, weightID)) %>%
#     mutate(lag = substr(parm, 4, (nchar(parm) - 1))) %>%
#     mutate(lag = as.numeric(lag)) %>%
#     mutate(variable = var) %>%
#     mutate(dataset = datasetID) %>%
#     rename(LCI_obs = `2.5%`,
#            median_obs = `50%`,
#            UCI_obs = `97.5%`) %>%
#     dplyr::select(lag, variable, dataset,
#                   LCI_obs, median_obs, UCI_obs) %>% 
#     mutate(medianobs_c = accumulate(median_obs, `+`)) %>%
#     rowwise() %>%
#     mutate(LCIobs_diff = medianobs_c - (median_obs - LCI_obs),
#            UCIobs_diff =  medianobs_c + (UCI_obs - median_obs)) %>%
#     mutate(UCIobs_diff = case_when(UCIobs_diff > 1 ~ 1,
#                                    TRUE ~ UCIobs_diff))
#   
#   df_all <- df %>%
#     left_join(df2, by = c("lag", "variable", 'dataset'))
#   
#   plot <- ggplot(df_all) +
#     geom_ribbon(aes(x = lag, ymin = LCIobs_diff,
#                     ymax = UCIobs_diff),
#                 fill = observed_col, alpha = 0.5) +
#     geom_line(aes(x = lag, y = medianobs_c), 
#               color = observed_col, alpha = 0.5) +
#     geom_pointrange(aes(x = lag, y = median_c,
#                         ymin = LCI_diff, 
#                         ymax = UCI_diff),
#                     color = modeled_col) +    
#     geom_line(aes(x = lag, y = median_c), 
#               color = modeled_col, alpha = 0.5,
#               linewidth = 0.6) +
#     labs(x = "Seasons into the past", 
#          y = "Cumulative seasonal weights \n posterior median and 95% BCI")
#   
#   return(plot)
#   
# }
# 
# 
# weights_fun2 <- function(model, weightID, nlag, var, modID){
#   df <- as.data.frame(model$quantiles) %>%
#     rownames_to_column(var = "parm") %>%
#     filter(str_detect(parm, weightID)) %>%
#     mutate(lag = substr(parm, 4, (nchar(parm) - 1))) %>%
#     mutate(lag = as.numeric(lag)) %>%
#     mutate(variable = var) %>%
#     mutate(model = modID) %>%
#     mutate(sig = case_when(`2.5%` > (1/nlag) ~ "sig",
#                            TRUE ~ "notsig"))
#   
#   return(df)
# }
# 
# #fish weights
# fish1 <- weights_fun2(model = fish_sam, 
#                       weightID = "wA", 
#                       nlag = 6,
#                       var = "Plant Biomass",
#                       modID = "fish")
# 
# fish2 <- weights_fun2(model = fish_sam, 
#                       weightID = 'wB', 
#                       nlag = 11,
#                       var = "Temperature",
#                        modID = "fish")
# 
# #birds
# bird1 <- weights_fun2(model = bird_sam, 
#                       weightID = "wA", 
#                       nlag = 6,
#                       var = "Temperature",
#                       modID = "bird")
# 
# bird2 <- weights_fun2(model = bird_sam, 
#                       weightID = "wB", 
#                       nlag = 6,
#                       var = "Precipitation", 
#                       modID = 'bird')
# 
# #hoppers
# hop1 <- weights_fun2(model = sev_sam,
#                      weightID = "wA",
#                      nlag = 6,
#                      var = "Temperature",
#                      modID = "hop")
# 
# hop2 <- weights_fun2(model = sev_sam,
#                      weightID = "wB",
#                      nlag = 6,
#                      var = "Precipitation",
#                      modID = 'hop')
# 
# hop3 <- weights_fun2(model = sev_sam,
#                      weightID = "wC",
#                      nlag = 11,
#                      var = "Plant Biomass",
#                      modID = 'hop')
# 
# #plants
# plant1 <- weights_fun2(model = plant_sam,
#                        weightID = "wA",
#                        nlag = 8,
#                        var = "VPD",
#                        modID = 'plant')
# 
# plant2 <- weights_fun2(model = plant_sam,
#                        weightID = "wB",
#                        nlag = 8,
#                        var = "Precipitation",
#                        modID = 'plant')
# 
# modeled_weights <- fish1 %>%
#   bind_rows(fish2, bird1, bird2, hop1, 
#             hop2, hop3, plant1, plant2)
# 
# plotA <- ggplot(modeled_weights) +
#   geom_tile(aes(x = lag, y = variable, fill = `50%`)) +
#   facet_grid(model~., scales = "free") +
#   scale_fill_viridis_b(limits = c(0, 0.8), 
#                        oob = scales::squish,
#                        direction = -1) +
#   labs(title = "Modeled data")
# 
# #fish weights
# fish1r <- weights_fun2(model = fish_sam_raw, 
#                       weightID = "wA", 
#                       nlag = 6,
#                       var = "Plant Biomass",
#                       modID = "fish")
# 
# fish2r <- weights_fun2(model = fish_sam_raw, 
#                       weightID = 'wB', 
#                       nlag = 11,
#                       var = "Temperature",
#                       modID = "fish")
# 
# #birds
# bird1r <- weights_fun2(model = bird_sam_raw, 
#                       weightID = "wA", 
#                       nlag = 6,
#                       var = "Temperature",
#                       modID = "bird")
# 
# bird2r <- weights_fun2(model = bird_sam_raw, 
#                       weightID = "wB", 
#                       nlag = 6,
#                       var = "Precipitation", 
#                       modID = 'bird')
# 
# #hoppers
# hop1r <- weights_fun2(model = sev_sam_raw,
#                      weightID = "wA",
#                      nlag = 6,
#                      var = "Temperature",
#                      modID = "hop")
# 
# hop2r <- weights_fun2(model = sev_sam_raw,
#                      weightID = "wB",
#                      nlag = 6,
#                      var = "Precipitation",
#                      modID = 'hop')
# 
# hop3r <- weights_fun2(model = sev_sam_raw,
#                      weightID = "wC",
#                      nlag = 11,
#                      var = "Plant Biomass",
#                      modID = 'hop')
# 
# #plants
# plant1r <- weights_fun2(model = plant_sam_raw,
#                        weightID = "wA",
#                        nlag = 8,
#                        var = "VPD",
#                        modID = 'plant')
# 
# plant2r <- weights_fun2(model = plant_sam_raw,
#                        weightID = "wB",
#                        nlag = 8,
#                        var = "Precipitation",
#                        modID = 'plant')
# 
# raw_weights <- fish1r %>%
#   bind_rows(fish2r, bird1r, bird2r, hop1r, 
#             hop2r, hop3r, plant1r, plant2r)
# 
# plotB <- ggplot(raw_weights) +
#   geom_tile(aes(x = lag, y = variable, fill = `50%`)) +
#   facet_grid(model~., scales = "free") +
#   scale_fill_viridis_b(limits = c(0, 0.8), 
#                        oob = scales::squish,
#                        direction = -1) +
#   labs(title = "Observed data")
# 
# plotB + plotA +
#   plot_layout(guides = "collect")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# weights_fun <- function(model_raw, model_model, weightID, nlags){
#   
#   df_raw <- as.data.frame(model_raw$quantiles) %>%
#     rownames_to_column(var = "parm") %>%
#     filter(str_detect(parm, weightID)) %>%
#     mutate(lag = substr(parm, 4, (nchar(parm) - 1))) %>%
#     mutate(lag = as.numeric(lag)) %>%
#     mutate(type = "observed")
#   
#   df_model <- as.data.frame(model_model$quantiles) %>%
#     rownames_to_column(var = "parm") %>%
#     filter(str_detect(parm, weightID)) %>%
#     mutate(lag = substr(parm, 4, (nchar(parm) - 1))) %>%
#     mutate(lag = as.numeric(lag)) %>%
#     mutate(type = "modeled")
#   
#   df <- df_raw %>% rbind(df_model)
#   
#   plot <- ggplot(df) +
#     geom_hline(yintercept = 1/nlags, linetype = 2, linewidth = 0.1) +
#     geom_pointrange(aes(x = lag, y = `50%`, 
#                         color = type, 
#                         ymin = `2.5%`, ymax = `97.5%`),
#                     position = position_dodge(width = 1), size = 0.4)
#   
#   return(plot)
#   
# }
# 
# #fish kelp
# fishkelp <- weights_fun(model_raw = fish_sam_raw,
#             model_model = fish_sam,
#             weightID = "wA",
#             nlags = 6)
# 
# #fish temp
# fishtemp <- weights_fun(model_raw = fish_sam_raw,
#             model_model = fish_sam,
#             weightID = "wB",
#             nlags = 11)
# 
# #birds
# #"Temperature", "Precipitation"
# #temp
# birdtemp <- weights_fun(model_raw = bird_sam_raw,
#             model_model = bird_sam,
#             weightID = "wA",
#             nlags = 6)
# #precip
# birdppt <- weights_fun(model_raw = bird_sam_raw,
#             model_model = bird_sam,
#             weightID = "wB",
#             nlags = 6)
# 
# #grasshoppers
# #"Temperature","Precipitation", "NPP"
# #temp
# hoptemp <- weights_fun(model_raw = sev_sam_raw,
#             model_model = sev_sam,
#             weightID = "wA",
#             nlags = 6)
# #ppt
# hopppt <- weights_fun(model_raw = sev_sam_raw,
#             model_model = sev_sam,
#             weightID = "wB",
#             nlags = 6)
# #npp
# hopnpp <- weights_fun(model_raw = sev_sam_raw,
#             model_model = sev_sam,
#             weightID = "wC",
#             nlags = 11)
# 
# #plants
# #"VPD","Precipitation"
# #vpd
# plantvpd <- weights_fun(model_raw = plant_sam_raw,
#             model_model = plant_sam,
#             weightID = "wA",
#             nlags = 8)
# #ppt
# plantppt <- weights_fun(model_raw = plant_sam_raw,
#             model_model = plant_sam,
#             weightID = "wB",
#             nlags = 8)
# 
# 
# # combine -----------------------------------------------------------------
# 
# fishkelp + fishtemp + birdtemp + birdppt +
#   hoptemp + hopppt + hopnpp + plantvpd + plantppt +
#   plot_layout(guides = "collect")
# 
# fishtemp + birdppt + hoptemp + hopppt + hopnpp+
#   plot_layout(guides = "collect")

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


# Wights plot -------------------------------------------------------------

# fish --------------------------------------------------------------------

#temperature
# fish_tweights_raw <- as.data.frame(fish_sam_raw$quantiles) %>%
#   rownames_to_column(var = "parm") %>%
#   filter(str_detect(parm, "wB")) %>%
#   mutate(season = case_when(parm %in% c("wB[1]", "wB[3]", "wB[5]",
#                                         'wB[7]', 'wB[9]','wB[11]') ~ "Warm",
#                             parm %in% c("wB[2]", "wB[4]", 'wB[6]',
#                                         'wB[8]', 'wB[10]') ~ "Cold")) %>%
#   mutate(year = case_when(parm == "wB[1]" ~ 0,
#                           parm %in% c("wB[2]", 'wB[3]') ~ 1,
#                           parm %in% c('wB[4]', 'wB[5]') ~ 2,
#                           parm %in% c("wB[6]", 'wB[7]') ~ 3,
#                           parm %in% c("wB[8]", 'wB[9]') ~ 4,
#                           parm %in% c("wB[10]", 'wB[11]') ~ 5,
#                           TRUE ~ NA_real_)) %>% 
#   mutate(season_num = case_when(parm == "wB[1]" ~ 1,
#                                 parm == "wB[2]" ~ 2, 
#                                 parm == 'wB[3]' ~ 3,
#                                 parm == 'wB[4]' ~ 4,
#                                 parm == 'wB[5]' ~ 5,
#                                 parm == "wB[6]" ~ 6,
#                                 parm == 'wB[7]' ~ 7,
#                                 parm == "wB[8]" ~ 8,
#                                 parm == 'wB[9]' ~ 9,
#                                 parm == "wB[10]" ~ 10,
#                                 parm == 'wB[11]' ~ 11,
#                                 TRUE ~ NA_real_)) %>% 
#   mutate(above = case_when(
#     `50%` > 1/11 ~ "yes",
#     TRUE ~ "no"
#   )) %>% 
#   complete(season, year) %>% 
#   unite("color", season, above, sep = "-", remove = FALSE) %>%
#   mutate(type = "raw")
# 
# fish_tweights <- as.data.frame(fish_sam$quantiles) %>%
#   rownames_to_column(var = "parm") %>%
#   filter(str_detect(parm, "wB")) %>%
#   mutate(season = case_when(parm %in% c("wB[1]", "wB[3]", "wB[5]",
#                                         'wB[7]', 'wB[9]','wB[11]') ~ "Warm",
#                             parm %in% c("wB[2]", "wB[4]", 'wB[6]',
#                                         'wB[8]', 'wB[10]') ~ "Cold")) %>%
#   mutate(year = case_when(parm == "wB[1]" ~ 0,
#                           parm %in% c("wB[2]", 'wB[3]') ~ 1,
#                           parm %in% c('wB[4]', 'wB[5]') ~ 2,
#                           parm %in% c("wB[6]", 'wB[7]') ~ 3,
#                           parm %in% c("wB[8]", 'wB[9]') ~ 4,
#                           parm %in% c("wB[10]", 'wB[11]') ~ 5,
#                           TRUE ~ NA_real_)) %>% 
#   mutate(season_num = case_when(parm == "wB[1]" ~ 1,
#                                 parm == "wB[2]" ~ 2, 
#                                 parm == 'wB[3]' ~ 3,
#                                 parm == 'wB[4]' ~ 4,
#                                 parm == 'wB[5]' ~ 5,
#                                 parm == "wB[6]" ~ 6,
#                                 parm == 'wB[7]' ~ 7,
#                                 parm == "wB[8]" ~ 8,
#                                 parm == 'wB[9]' ~ 9,
#                                 parm == "wB[10]" ~ 10,
#                                 parm == 'wB[11]' ~ 11,
#                                 TRUE ~ NA_real_)) %>% 
#   mutate(above = case_when(
#     `50%` > 1/11 ~ "yes",
#     TRUE ~ "no"
#   )) %>% 
#   complete(season, year) %>% 
#   unite("color", season, above, sep = "-", remove = FALSE)%>%
#   mutate(type = "modeled")
# 
# all_fish_twt <- fish_tweights %>%
#   rbind(fish_tweights_raw)
# 
# ggplot(all_fish_twt) +
#   geom_rect(aes(xmin = 1.5, xmax = 2.5, ymin = 0, ymax = 0.55), 
#             fill = "gray92") +
#   geom_rect(aes(xmin = 3.5, xmax = 4.5, ymin = 0, ymax = 0.55), 
#             fill = "gray92") +
#   geom_rect(aes(xmin = 5.5, xmax = 6.5, ymin = 0, ymax = 0.55), 
#             fill = "gray92") +
#   geom_rect(aes(xmin = 7.5, xmax = 8.5, ymin = 0, ymax = 0.55), 
#             fill = "gray92") +
#   geom_rect(aes(xmin = 9.5, xmax = 10.5, ymin = 0, ymax = 0.55), 
#             fill = "gray92") +
#   geom_hline(yintercept = 1/11, linetype = 2, linewidth = 0.1) +
#   geom_pointrange(aes(x = season_num, y = `50%`, 
#                       color = type, 
#                       ymin = `2.5%`, ymax = `97.5%`),
#                   position = position_dodge(width = 1), size = 0.4) +
#   scale_x_continuous(breaks = c(1, 2, 4,  6,  8,  10), 
#                      labels = c(0, 1, 2,  3,  4, 5 )) +
#   scale_color_manual(values = c(modeled_col, observed_col)) +
#   labs(x = "Years into past",
#        y = "Importance weight \n (median and 95% BCI)")

