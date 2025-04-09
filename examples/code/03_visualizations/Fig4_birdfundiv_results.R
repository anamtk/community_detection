#bird fundiv SAM model results
#Ana Miller-ter Kuile
#February 6, 2025

#summariseing results from bird functional richness models

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

#Q
q_raw <- readRDS(here('examples',
                      'data_output',
                      'bird_fundiv',
                      'model_outputs',
                      'bird_q_betareg_summary_raw.RDS'))

q_corr <- readRDS(here('examples',
                      'data_output',
                      'bird_fundiv',
                      'model_outputs',
                      'bird_q_betareg_summary_corrected.RDS'))

data_corrected <- readRDS(here('examples',
                    "data_output",
                    'bird_fundiv',
                    'tidy_data',
                    'bird_fd_metrics_corrected.RDS'))

data_raw <- readRDS(here('examples',
                               "data_output",
                               'bird_fundiv',
                               'tidy_data',
                               'bird_fd_metrics_raw.RDS'))

data_all <- data_corrected %>%
  left_join(data_raw, by =c("TransID", "yrID"))

# Raw corrected values ----------------------------------------------------
#how do raw/correcetd compare?
data_rao <- data_all %>%
  dplyr::select(TransID, yrID, Q_mean, Q) %>%
  rename(modeled = Q_mean,
         observed = Q) %>%
  pivot_longer(modeled:observed,
               names_to = "type",
               values_to = "Q")

(raodat_plot <- ggplot(data_rao, aes(x = type,
                     fill = type,
                     y = Q)) +
  geom_boxplot()+
  scale_fill_manual(values = c(modeled_col, observed_col),
                     labels = c("Modeled", "Empirical")) +
  theme(legend.position = "none") +
    labs(x = "Data type",
         y = "Rao's quadratic entropy"))

# sam plot function -------------------------------------------------------

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

# (rao_plots <- sam_plot_fun(model_model = q_corr,
#                           model_raw = q_raw)+
#     scale_y_discrete(labels = c("Temperature", "Precipitation")) + 
#     labs(title = "Bird Rao"))

rao_df <- sam_sum_fun(model_model = q_corr,
             model_raw = q_raw)

(modeled_raoeff <- ggplot(rao_df) +
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

(observed_raoeff <- ggplot(rao_df) +
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
    scale_color_manual(values = c(modeled_light, observed_col),
                       labels = c("Modeled", "Empirical")) +
    scale_shape_manual(values = c('nonsig' = 21, 'sig' = 19)) +
    scale_y_discrete(labels = c("Temperature", "Precipitation")) + 
    labs(x = "Covariate effect \n (Median and 95% BCI)", y = "") +
    theme(legend.position = "none"))

(raoeff <- ggplot(rao_df) +
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
    scale_color_manual(values = c(modeled_col, observed_col),
                       labels = c("Modeled", "Empirical")) +
    scale_shape_manual(values = c('nonsig' = 21, 'sig' = 19)) +
    scale_y_discrete(labels = c("Temperature", "Precipitation")) + 
    labs(x = "Covariate effect \n (Median and 95% BCI)", y = "") +
    theme(legend.position = "none"))

# Rao weight plots ---------------------------------------------------------------
rao_raw_w <- as.data.frame(q_raw$quantiles) %>%
  rownames_to_column(var = "par") %>%
  filter(str_detect(par, "w")) %>%
  filter(!str_detect(par, "cumm")) %>%
  mutate(lag = str_sub(par, 4, 4)) %>%
  mutate(covariate = case_when(str_detect(par, "wA") ~ "Temperature",
                               str_detect(par, "wB") ~ "Precipitation")) %>%
  dplyr::select(lag, covariate, `2.5%`, `50%`, `97.5%`) %>%
  mutate(type = "observed")

rao_corr_w <- as.data.frame(q_corr$quantiles) %>%
  rownames_to_column(var = "par") %>%
  filter(str_detect(par, "w")) %>%
  filter(!str_detect(par, "cumm")) %>%
  mutate(lag = str_sub(par, 4, 4)) %>%
  mutate(covariate = case_when(str_detect(par, "wA") ~ "Temperature",
                               str_detect(par, "wB") ~ "Precipitation")) %>%
  dplyr::select(lag, covariate, `2.5%`, `50%`, `97.5%`) %>%
  mutate(type = "modeled")

rao_ws <- rao_raw_w %>%
  bind_rows(rao_corr_w)

(rao_weight_plot <- ggplot(rao_ws, aes(x = lag, y = `50%`, color = type)) +
    geom_hline(yintercept = 1/4, linetype = 2, alpha = 0.4) +
    geom_pointrange(aes(x = lag, 
                        y = `50%`,
                        ymin = `2.5%`,
                        ymax = `97.5%`,
                        color = type),
                    position=position_dodge(width=0.5),
                    size = 0.4) +
    facet_grid(~covariate)+
    scale_color_manual(values = c(modeled_col, observed_col),
                       labels = c("Modeled", "Empirical")) +
    labs(x= "Seasons into past",
         y = "Importance weight\n(median and 95% BCI)"))

rao_raw_w2 <- as.data.frame(q_raw$quantiles) %>%
  rownames_to_column(var = "par") %>%
  filter(str_detect(par, "cumm")) %>%
  mutate(lag = str_sub(par, (nchar(par)-1), -2)) %>%
  mutate(covariate = case_when(str_detect(par, "temp") ~ "Temperature",
                               str_detect(par, "ppt") ~ "Precipitation")) %>%
  dplyr::select(lag, covariate, `2.5%`, `50%`, `97.5%`) %>%
  rename(LCI_obs = `2.5%`,
         median_obs = `50%`,
         UCI_obs = `97.5%`) 

rao_corr_w2 <- as.data.frame(q_corr$quantiles) %>%
  rownames_to_column(var = "par") %>%
  filter(str_detect(par, "cumm")) %>%
  mutate(lag = str_sub(par, (nchar(par)-1), -2)) %>%
  mutate(covariate = case_when(str_detect(par, "temp") ~ "Temperature",
                               str_detect(par, 'ppt') ~ "Precipitation")) %>%
  dplyr::select(lag, covariate, `2.5%`, `50%`, `97.5%`) 

rao_ws2 <- rao_raw_w2 %>%
  left_join(rao_corr_w2, by = c("lag", "covariate")) %>%
  mutate(lag = as.numeric(lag))

rao_weights_plot2 <- rao_ws2 %>%
  ggplot() +
  geom_ribbon(aes(x = lag, 
                  ymin = LCI_obs,
                  ymax = UCI_obs),
              fill = observed_col, 
              alpha = 0.5) +
  geom_line(aes(x = lag, 
                y = median_obs), 
            color = observed_col,
            alpha = 0.5) +
  geom_pointrange(aes(x = lag, 
                      y = `50%`,
                      ymin = `2.5%`,
                      ymax = `97.5%`),
                  color = modeled_col,
                  size = 0.4) +  
  geom_line(aes(x = lag, 
                y = `50%`), 
            color = modeled_col,
            alpha = 0.5,
            linewidth = 0.6) +
  facet_wrap(~covariate) +
  labs(x = "Seasons into past",
       y = "Cumulative seasonal weights \n posterior median and 95% BCI")

(rao_weights_precip_mod <- rao_ws2 %>%
    filter(covariate == "Precipitation") %>%
  ggplot() +
  geom_ribbon(aes(x = lag, 
                  ymin = LCI_obs,
                  ymax = UCI_obs),
              fill = observed_col, 
              alpha = 0.2) +
  geom_line(aes(x = lag, 
                y = median_obs), 
            color = observed_col,
            alpha = 0.2) +
  geom_pointrange(aes(x = lag, 
                      y = `50%`,
                      ymin = `2.5%`,
                      ymax = `97.5%`),
                  color = modeled_col,
                  size = 0.4) +  
  geom_line(aes(x = lag, 
                y = `50%`), 
            color = modeled_col,
            alpha = 0.5,
            linewidth = 0.6) +
  labs(x = "Seasons into past",
       y = "Cumulative seasonal weights \n posterior median and 95% BCI",
       title = "Precipitation")
)

(rao_weights_precip_obs <- rao_ws2 %>%
    filter(covariate == "Precipitation") %>%
    ggplot() +
    geom_ribbon(aes(x = lag, 
                    ymin = `2.5%`,
                    ymax = `97.5%`),
                fill = modeled_col, 
                alpha = 0.2) +
    geom_line(aes(x = lag, 
                  y = `50%`), 
              color = modeled_col,
              alpha = 0.2) +
    geom_pointrange(aes(x = lag, 
                        y = median_obs,
                        ymin = LCI_obs,
                        ymax = UCI_obs),
                    color = observed_col,
                    size = 0.4) +  
    geom_line(aes(x = lag, 
                  y = median_obs), 
              color = observed_col,
              alpha = 0.5,
              linewidth = 0.6) +
    labs(x = "Seasons into past",
         y = "Cumulative seasonal weights \n posterior median and 95% BCI",
         title = "Precipitation")
)

(rao_weights_temp_mod <- rao_ws2 %>%
    filter(covariate == "Temperature") %>%
    ggplot() +
    geom_pointrange(aes(x = lag, 
                        y = `50%`,
                        ymin = `2.5%`,
                        ymax = `97.5%`),
                    color = modeled_col,
                    size = 0.4) +  
    geom_line(aes(x = lag, 
                  y = `50%`), 
              color = modeled_col,
              alpha = 0.5,
              linewidth = 0.6) +
    labs(x = "Seasons into past",
         y = "",
         title = "Temperature")
)
# Plot --------------------------------------------------------------------


rao_weights_precip_mod +rao_weights_temp_mod +
rao_weights_precip_obs

ggsave(here('examples',
            'pictures',
            'original_R',
            'Fig3_weights_plots.pdf'),
       height =3,
       width = 7.5,
       units = "in")
raodat_plot


ggsave(here('examples',
            'pictures',
            'original_R',
            'Fig3_data_boxplot.pdf'),
       height =3,
       width = 3,
       units = "in")
raoeff

ggsave(here('examples',
            'pictures',
            'original_R',
            'Fig3_effects_plot.pdf'),
       height =3,
       width = 4,
       units = "in")

observed_raoeff
ggsave(here('examples',
            'pictures',
            'original_R',
            'Fig3_observedeffects_plot.pdf'),
       height =2,
       width = 3,
       units = "in")

modeled_raoeff
ggsave(here('examples',
            'pictures',
            'original_R',
            'Fig3_modeledeffects_plot.pdf'),
       height =2,
       width = 3,
       units = "in")

rao_weights_precip_mod +rao_weights_temp_mod 

ggsave(here('examples',
            'pictures',
            'original_R',
            'Fig3_modeled_weights_plots.pdf'),
       height =3,
       width = 5,
       units = "in")

rao_weights_precip_obs

ggsave(here('examples',
            'pictures',
            'original_R',
            'Fig3_observed_weights_plots.pdf'),
       height =3,
       width = 3,
       units = "in")
