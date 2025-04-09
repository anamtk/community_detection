#plant gain/loss SAM model results
#Ana Miller-ter Kuile
#February 6, 2025

#summariseing results from plant stability models


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

source(here('examples',
            'code',
            '00_functions',
            'tidy_functions.R'))

source(here('examples',
            'code',
            '00_functions',
            'plot_functions.R'))

modeled_col <- "#E88C23"
observed_col <- "#438AA8"

observed_light <- "#A9C7D5"
modeled_light <- "#F4C79F"

# Load data ---------------------------------------------------------------
#BIRD DATASETS
#FRIC
fric_raw <- readRDS(here('examples',
                         'data_output',
                         'bird_fundiv',
                         'model_outputs',
                         'bird_fric_betareg_summary_raw.RDS'))

fric_raw_cumul <- readRDS(here('examples',
                               'data_output',
                               'bird_fundiv',
                               'model_outputs',
                               'bird_fric_betareg_summary_raw_cumulativewts.RDS'))

fric_corr <- readRDS(here('examples',
                          'data_output',
                          'bird_fundiv',
                          'model_outputs',
                          'bird_fric_betareg_summary_corrected.RDS'))

bird_corrected <- readRDS(here('examples',
                               "data_output",
                               'bird_fundiv',
                               'tidy_data',
                               'bird_fd_metrics_corrected.RDS'))

bird_raw <- readRDS(here('examples',
                         "data_output",
                         'bird_fundiv',
                         'tidy_data',
                         'bird_fd_metrics_raw.RDS'))

bird_all <- bird_corrected %>%
  left_join(bird_raw, by =c("TransID", "yrID"))

#PLANT DATASETS
gain_raw <- readRDS(here('examples',
                         "data_output",
                         'plant_turnover',
                         'model_outputs',
                         'plant_gains_betareg_summary_raw.RDS'))

gain_corr <- readRDS(here('examples',
                          "data_output",
                          'plant_turnover',
                          'model_outputs',
                          'plant_gains_betareg_summary_corrected.RDS'))

plant_raw <- read.csv(here('examples',
                          'data_output',
                          'plant_turnover',
                          'tidy_data',
                          'plant_betareg_tidydata.csv'))


# raw_corrected comparison ------------------------------------------------

#how do raw/correcetd compare?
data_fric <- bird_all %>%
  dplyr::select(TransID, yrID, FRic_mean, FRic) %>%
  rename(modeled = FRic_mean,
         observed = FRic) %>%
  pivot_longer(modeled:observed,
               names_to = "type",
               values_to = "FRic")

(fricdat_plot <- ggplot(data_fric, aes(x = type, 
                                       fill = type,
                                       y = FRic)) +
    geom_boxplot()+
    scale_fill_manual(values = c(modeled_col, observed_col),
                      labels = c("Modeled", "Empirical"))+
    theme(legend.position = "none")) 

data_gain <- plant_raw %>%
  dplyr::select(gain, mean_gain) %>%
  rename(modeled = mean_gain,
         observed = gain) %>%
  pivot_longer(modeled:observed,
               names_to = "type",
               values_to = "gain")

(gain_plot_box <- ggplot(data_gain, aes(x = type, y = gain, fill = type)) +
    geom_boxplot()+
    scale_fill_manual(values = c(modeled_col, observed_col),
                      labels = c("Modeled", "Empirical"))+
    labs(x = "Data type",
         y = "Species gains") +
    theme(legend.position = "none"))



# SAM results -------------------------------------------------------------

sam_sum_fun <- function(model_model, model_raw){
  
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
                    show.legend = T,
                    position=position_dodge(width=0.5),
                    size = 0.4) +
    guides(shape = FALSE) +
    scale_color_manual(values = c(modeled_col, observed_col),
                       labels = c("Modeled", "Empirical")) +
    scale_shape_manual(values = c('nonsig' = 1, 'sig' = 19)) +
    labs(x = "Covariate effect \n (Median and 95% BCI)", y = "") 
  
  return(all_beta)
  
}


weights_fun <- function(model_model, model_raw){
  raw_w <- as.data.frame(model_raw$quantiles) %>%
    rownames_to_column(var = "par") %>%
    filter(str_detect(par, "w")) %>%
    mutate(lag = str_sub(par, 4, (nchar(par)-1))) %>%
    mutate(covariate = case_when(str_detect(par, "wA") ~ "VPD",
                                 str_detect(par, "wB") ~ "Precipitation")) %>%
    dplyr::select(lag, covariate, `2.5%`, `50%`, `97.5%`) %>%
    mutate(type = "observed")
  
  corr_w <- as.data.frame(model_model$quantiles) %>%
    rownames_to_column(var = "par") %>%
    filter(str_detect(par, "w")) %>%
    filter(!str_detect(par, "cumm")) %>%
    mutate(lag = str_sub(par, 4, 4)) %>%
    mutate(covariate = case_when(str_detect(par, "wA") ~ "VPD",
                                 str_detect(par, "wB") ~ "Precipitation")) %>%
    dplyr::select(lag, covariate, `2.5%`, `50%`, `97.5%`) %>%
    mutate(type = "modeled")
  
  ws <- raw_w %>%
    bind_rows(corr_w)
  
  weight_plot <- ggplot(ws, aes(x = lag, y = `50%`, color = type)) +
    geom_hline(yintercept = 1/4, linetype = 2, alpha = 0.4) +
    geom_pointrange(aes(x = lag, 
                        y = `50%`,
                        ymin = `2.5%`,
                        ymax = `97.5%`,
                        color = type),
                    position=position_dodge(width=0.5),
                    size = 0.4) +
    facet_grid(covariate~.)+
    scale_color_manual(values = c(modeled_col, observed_col),
                       labels = c("Modeled", "Empirical")) +
    labs(x= "Seasons into past",
         y = "Importance weight\n(median and 95% BCI)")
  
  return(weight_plot)
  
}


# Bird FRIC ---------------------------------------------------------------

bird_df<- sam_sum_fun(fric_corr, fric_raw)

(modeled_bird <- ggplot(bird_df) +
    geom_vline(xintercept = 0, linetype = 2, alpha = 0.4) +
    geom_errorbar(aes(xmin = `2.5%`,
                      xmax = `97.5%`,
                      y = parm,
                      color = type),
                  width = 0,
                  position=position_dodge(width=0.5)) +
    geom_point(aes(x = `50%`,
                   y = parm, 
                   color = type, 
                   fill = type,
                   shape = sig), 
               show.legend = T,
               position=position_dodge(width=0.5),
               size = 4) +
    guides(shape = FALSE,
           fill = FALSE) +
    scale_fill_manual(values = c('white', 'white'),
                      labels = c("Modeled", "Empirical")) +
    scale_color_manual(values = c(modeled_col, observed_light),
                       labels = c("Modeled", "Empirical")) +
    scale_shape_manual(values = c('nonsig' = 21, 'sig' = 19)) +
    scale_y_discrete(labels = c("Temperature", "Precipitation")) + 
    labs(x = "Covariate effect \n (Median and 95% BCI)", y = "") +
    theme(legend.position = "none"))

corr_temp <- as.data.frame(fric_corr$quantiles)%>%
  rownames_to_column(var = "par") %>%
  filter(str_detect(par, "wA")) %>%
  filter(!str_detect(par, "cumm")) %>%
  mutate(lag = str_sub(par, 4, 4)) %>%
  mutate(covariate = case_when(str_detect(par, "wA") ~ "Temperature")) %>%
  dplyr::select(lag, covariate, `2.5%`, `50%`, `97.5%`) %>%
  mutate(type = "modeled")

(temp_plot <- ggplot(corr_temp, aes(x = lag, y = `50%`)) +
    geom_hline(yintercept = 1/8, linetype = 2, alpha = 0.4) +
    geom_pointrange(aes(x = lag, 
                        y = `50%`,
                        ymin = `2.5%`,
                        ymax = `97.5%`),
                    color = modeled_col,
                    position=position_dodge(width=0.5),
                    size = 0.4) +
    facet_grid(covariate~.)+
    labs(x= "Seasons into past",
         y = "Importance weight\n(median and 95% BCI)"))

raw_ppt <- as.data.frame(fric_raw$quantiles)%>%
  rownames_to_column(var = "par") %>%
  filter(str_detect(par, "wB")) %>%
  filter(!str_detect(par, "cumm")) %>%
  mutate(lag = str_sub(par, 4, 4)) %>%
  mutate(covariate = case_when(str_detect(par, "wB") ~ "Precipitation")) %>%
  dplyr::select(lag, covariate, `2.5%`, `50%`, `97.5%`) %>%
  mutate(type = "observed")


(ppt_plot <- ggplot(raw_ppt, aes(x = lag, y = `50%`)) +
    geom_hline(yintercept = 1/8, linetype = 2, alpha = 0.4) +
    geom_pointrange(aes(x = lag, 
                        y = `50%`,
                        ymin = `2.5%`,
                        ymax = `97.5%`),
                    color = observed_light,
                    position=position_dodge(width=0.5),
                    size = 0.4) +
    facet_grid(covariate~.)+
    labs(x= "Seasons into past",
         y = "Importance weight\n(median and 95% BCI)"))

modeled_bird + (temp_plot / ppt_plot) +
  plot_annotation(tag_levels = "A")

# Plant gains -------------------------------------------------------------


gain_df <- sam_sum_fun(model_model = gain_corr,
                         model_raw = gain_raw)

(modeled_gaineff <- ggplot(gain_df) +
    geom_vline(xintercept = 0, linetype = 2, alpha = 0.4) +
    geom_errorbar(aes(xmin = `2.5%`,
                      xmax = `97.5%`,
                      y = parm,
                      color = type),
                  width = 0,
                  position=position_dodge(width=0.5)) +
    geom_point(aes(x = `50%`,
                   y = parm, 
                   color = type, 
                   fill = type,
                   shape = sig), 
               show.legend = T,
               position=position_dodge(width=0.5),
               size = 4) +
    guides(shape = FALSE,
           fill = FALSE) +
    scale_fill_manual(values = c('white', 'white'),
                      labels = c("Modeled", "Empirical")) +
    scale_color_manual(values = c(modeled_col, observed_light),
                       labels = c("Modeled", "Empirical")) +
    scale_shape_manual(values = c('nonsig' = 21, 'sig' = 19)) +
    scale_y_discrete(labels = c("VPD", "Precipitation")) + 
    labs(x = "Covariate effect \n (Median and 95% BCI)", y = "") +
    theme(legend.position = "none"))

corr_w <- as.data.frame(gain_corr$quantiles) %>%
  rownames_to_column(var = "par") %>%
  filter(str_detect(par, "w")) %>%
  filter(!str_detect(par, "cumm")) %>%
  mutate(lag = str_sub(par, 4, 4)) %>%
  mutate(covariate = case_when(str_detect(par, "wA") ~ "VPD",
                               str_detect(par, "wB") ~ "Precipitation")) %>%
  dplyr::select(lag, covariate, `2.5%`, `50%`, `97.5%`) %>%
  mutate(type = "modeled")

(weight_plot <- ggplot(corr_w, aes(x = lag, y = `50%`)) +
    geom_hline(yintercept = 1/8, linetype = 2, alpha = 0.4) +
    geom_pointrange(aes(x = lag, 
                        y = `50%`,
                        ymin = `2.5%`,
                        ymax = `97.5%`),
                    color = modeled_col,
                    position=position_dodge(width=0.5),
                    size = 0.4) +
    facet_grid(covariate~.)+
    labs(x= "Seasons into past",
         y = "Importance weight\n(median and 95% BCI)"))


# Export ------------------------------------------------------------------

fricdat_plot + gain_plot_box +
  plot_annotation(tag_levels = "A")

ggsave(here('examples',
            'pictures',
            'original_R',
            'SI_other_metrics_boxplots.pdf'),
       height =3,
       width = 5,
       units = "in")

modeled_bird + (temp_plot / ppt_plot) +
  plot_annotation(tag_levels = "A")

ggsave(here('examples',
            'pictures',
            'original_R',
            'SI_bird_FRic.jpg'),
       height =4,
       width = 6,
       units = "in")

modeled_gaineff + weight_plot+
  plot_annotation(tag_levels = "A")

ggsave(here('examples',
            'pictures',
            'original_R',
            'SI_plant_gains.jpg'),
       height =4,
       width = 6,
       units = "in")
