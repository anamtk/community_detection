#grasshopper stability SAM model results
#Ana Miller-ter Kuile
#February 6, 2025

#summariseing results from grasshopper models

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

hop_raw <- readRDS(here('examples',
                        'data_output',
                        'grasshopper_stability',
                        'computing_cluster_outputs',
                        'sev_SAM_raw_summary.RDS'))

hop_raw_ppt <- readRDS(here('examples',
                            'data_output',
                            'grasshopper_stability',
                            'computing_cluster_outputs',
                            'sev_SAM_cummppt_raw.RDS')) 

hop_raw_temp <- readRDS(here('examples',
                            'data_output',
                            'grasshopper_stability',
                            'computing_cluster_outputs',
                            'sev_SAM_cummtemp_raw.RDS')) 

hop_raw_npp <- readRDS(here('examples',
                             'data_output',
                             'grasshopper_stability',
                             'computing_cluster_outputs',
                             'sev_SAM_cummnpp_raw.RDS')) 

hop_corr <- readRDS(here('examples',
                         'data_output',
                         'grasshopper_stability',
                         'computing_cluster_outputs',
                         'sev_SAM_summary.RDS'))

hop_corr_ppt <- readRDS(here('examples',
                            'data_output',
                            'grasshopper_stability',
                            'computing_cluster_outputs',
                            'sev_SAM_cummppt.RDS'))

hop_corr_temp <- readRDS(here('examples',
                             'data_output',
                             'grasshopper_stability',
                             'computing_cluster_outputs',
                             'sev_SAM_cummtemp.RDS')) 

hop_corr_npp <- readRDS(here('examples',
                            'data_output',
                            'grasshopper_stability',
                            'computing_cluster_outputs',
                            'sev_SAM_cummnpp.RDS'))

#also have cumulative summaries for both if we want that for graphing

data<- read.csv(here('examples',
                         'data_output',
                         'grasshopper_stability',
                         'tidy_data',
                         "grasshopper_betareg_tidydata.csv"))

# Data --------------------------------------------------------------------

data2 <- data %>%
  dplyr::select(mean_bray, observed_all) %>%
  rename('modeled' = mean_bray,
         observed = observed_all) %>%
  pivot_longer(modeled:observed,
               names_to = "type",
               values_to = 'bray')

(braydat_plot <- ggplot(data2, aes(x = type,
                                     fill = type,
                                     y = bray)) +
    geom_boxplot()+
    scale_fill_manual(values = c(modeled_col, observed_col),
                      labels = c("Modeled", "Empirical")) +
    theme(legend.position = "none") +
    labs(x = "Data type",
         y = "Bray-curtis temporal stability"))


# Sam betas ---------------------------------------------------------------


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

hop_betas <- sam_sum_fun(model_model = hop_corr,
                         model_raw = hop_raw)

(modeled_brayeff <- ggplot(hop_betas) +
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
    scale_y_discrete(labels = c("Temperature", "Precipitation", "Plant biomass")) + 
    labs(x = "Covariate effect \n (Median and 95% BCI)", y = "") +
    theme(legend.position = "none"))

# weights -----------------------------------------------------------------

hop_raw_w <- as.data.frame(hop_raw$quantiles) %>%
  rownames_to_column(var = "par") %>%
  filter(str_detect(par, "w")) %>%
  filter(!str_detect(par, "cumm")) %>%
  filter(!str_detect(par, "b0.web")) %>%
  filter(!str_detect(par, "sig.web")) %>%
  mutate(lag = str_sub(par, 4, (nchar(par)-1))) %>%
  mutate(covariate = case_when(str_detect(par, "wA") ~ "Temperature",
                               str_detect(par, "wB") ~ "Precipitation",
                               str_detect(par, "wC") ~ "Plant biomass")) %>%
  dplyr::select(lag, covariate, `2.5%`, `50%`, `97.5%`) %>%
  mutate(type = "observed")

hop_corr_w <- as.data.frame(hop_corr$quantiles) %>%
  rownames_to_column(var = "par") %>%
  filter(str_detect(par, "w")) %>%
  filter(!str_detect(par, "cumm")) %>%
  filter(!str_detect(par, "b0.web")) %>%
  filter(!str_detect(par, "sig.web")) %>%
  mutate(lag = str_sub(par, 4, (nchar(par)-1))) %>%
  mutate(covariate = case_when(str_detect(par, "wA") ~ "Temperature",
                               str_detect(par, "wB") ~ "Precipitation",
                               str_detect(par, "wC") ~ "Plant biomass")) %>%
  dplyr::select(lag, covariate, `2.5%`, `50%`, `97.5%`) %>%
  mutate(type = "modeled")

hop_ws <- hop_raw_w %>%
  bind_rows(hop_corr_w)


(hop_weight_plot1 <- hop_ws %>%
    filter(covariate == "Plant biomass") %>%
    mutate(lag = factor(lag, levels = c("1", "2", 
                                        "3", "4", "5",
                                        "6", "7", "8",
                                        "9", "10",
                                        "11", "12"))) %>%
    ggplot(aes(x = lag, y = `50%`, color = type)) +
    geom_hline(yintercept = 1/11, linetype = 2, alpha = 0.4) +
    geom_pointrange(aes(x = lag, 
                        y = `50%`,
                        ymin = `2.5%`,
                        ymax = `97.5%`,
                        color = type),
                    position=position_dodge(width=0.5),
                    size = 0.4) +
    facet_grid(~covariate)+
    scale_color_manual(values = c(modeled_col, observed_light ),
                       labels = c("Modeled", "Empirical")) +
    labs(x= "Seasons into past",
         y = "Importance weight\n(median and 95% BCI)"))

(hop_weight_plot2 <- hop_ws %>%
    filter(covariate != "Plant biomass") %>%
    mutate(lag = factor(lag, levels = c("1", "2", 
                                        "3", "4", "5",
                                        "6", "7", "8"))) %>%
    ggplot(aes(x = lag, y = `50%`, color = type)) +
    geom_hline(yintercept = 1/6, linetype = 2, alpha = 0.4) +
    #geom_line(aes(x = lag, y = `50%`, group = type, color = type)) +
    geom_pointrange(aes(x = lag, 
                        y = `50%`,
                        ymin = `2.5%`,
                        ymax = `97.5%`,
                        color = type),
                    position=position_dodge(width=0.5),
                    size = 0.4) +
    facet_grid(~covariate)+
    scale_color_manual(values = c(modeled_col, observed_light ),
                       labels = c("Modeled", "Empirical")) +
    labs(x= "Seasons into past",
         y = "Importance weight\n(median and 95% BCI)")+
    theme(legend.position = "none"))

(hop_weight_plot3 <- hop_ws %>%
    filter(covariate != "Plant biomass") %>%
    mutate(lag = factor(lag, levels = c("1", "2", 
                                        "3", "4", "5",
                                        "6", "7", "8"))) %>%
    ggplot(aes(x = lag, y = `50%`, color = type)) +
    geom_hline(yintercept = 1/8, linetype = 2, alpha = 0.4) +
    geom_ribbon(aes(x = lag, 
                    ymin = `2.5%`, 
                    ymax = `97.5%`,
                    group = type,
                    fill = type),
                alpha = 0.2,
                color = "white") +
    geom_line(aes(x = lag, y = `50%`, group = type, color = type)) +
    # geom_pointrange(aes(x = lag, 
    #                     y = `50%`,
    #                     ymin = `2.5%`,
    #                     ymax = `97.5%`,
    #                     color = type),
    #                 position=position_dodge(width=0.5),
    #                 size = 0.4) +
    facet_grid(~covariate)+
    scale_color_manual(values = c(modeled_col, observed_col),
                       labels = c("Modeled", "Empirical")) +
    scale_fill_manual(values = c(modeled_col, observed_col),
                       labels = c("Modeled", "Empirical")) +
    labs(x= "Seasons into past",
         y = "Importance weight\n(median and 95% BCI)"))


hop_ppt_wt <- hop_raw_ppt %>%
  left_join(hop_corr_ppt, by = c("lag", "variable", "dataset"))

ggplot(hop_ppt_wt) +
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
                      y = median,
                      ymin = LCI,
                      ymax = UCI),
                  color = modeled_col,
                  size = 0.4) +  
  geom_line(aes(x = lag, 
                y = median), 
            color = modeled_col,
            alpha = 0.5,
            linewidth = 0.6) +
  labs(x = "Seasons into past",
       y = "Cumulative seasonal weights \n posterior median and 95% BCI",
       title = "Precipitation")

hop_temp_wt <- hop_raw_temp %>%
  left_join(hop_corr_temp, by = c("lag", "variable", "dataset"))

ggplot(hop_temp_wt) +
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
                      y = median,
                      ymin = LCI,
                      ymax = UCI),
                  color = modeled_col,
                  size = 0.4) +  
  geom_line(aes(x = lag, 
                y = median), 
            color = modeled_col,
            alpha = 0.5,
            linewidth = 0.6) +
  labs(x = "Seasons into past",
       y = "Cumulative seasonal weights \n posterior median and 95% BCI",
       title = "Temperature")

hop_npp_wt <- hop_raw_npp %>%
  left_join(hop_corr_npp, by = c("lag", "variable", "dataset"))

(hopnppwtplot <- hop_npp_wt %>%
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
                      y = median,
                      ymin = LCI,
                      ymax = UCI),
                  color = modeled_col,
                  size = 0.4) +  
  geom_line(aes(x = lag, 
                y = median), 
            color = modeled_col,
            alpha = 0.5,
            linewidth = 0.6) +
  labs(x = "Seasons into past",
       y = "Cumulative seasonal weights \n posterior median and 95% BCI",
       title = "Plant biomass") +
    scale_x_continuous(breaks = c(1:11)))


# Together ----------------------------------------------------------------

modeled_brayeff

ggsave(here('examples',
            'pictures',
            'original_R',
            'Fig5_grasshoppermodeledeffects_plot.pdf'),
       height =2,
       width = 3,
       units = "in")

hopnppwtplot /
hop_weight_plot2
ggsave(here('examples',
            'pictures',
            'original_R',
            'Fig5_modeled_weights_plots.pdf'),
       height =5,
       width = 4,
       units = "in")

braydat_plot 
