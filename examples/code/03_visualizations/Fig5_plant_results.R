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

loss_raw <- readRDS(here('examples',
                         "data_output",
                         'plant_turnover',
                         'model_outputs',
                         'plant_loss_betareg_summary_raw.RDS'))

loss_corr <- readRDS(here('examples',
                          "data_output",
                          'plant_turnover',
                          'model_outputs',
                          'plant_loss_betareg_summary_corrected.RDS'))

data_raw <- read.csv(here('examples',
                          'data_output',
                          'plant_turnover',
                          'tidy_data',
                          'plant_betareg_tidydata.csv'))


# raw_corrected comparison ------------------------------------------------

data_loss <- data_raw %>%
  dplyr::select(loss, mean_loss) %>%
  rename(modeled = mean_loss,
         observed = loss) %>%
  pivot_longer(modeled:observed,
               names_to = "type",
               values_to = "loss")

(loss_plot_box <- ggplot(data_loss, aes(x = type, y = loss, fill = type)) +
  geom_boxplot()+
  scale_fill_manual(values = c(modeled_col, observed_col),
                    labels = c("Modeled", "Empirical"))+
    labs(x = "Data type",
         y = "Species losses") +
    theme(legend.position = "none"))
  
# Effects plots -----------------------------------------------------------

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


loss_df <- sam_sum_fun(model_model = loss_corr,
                      model_raw = loss_raw)

(modeled_losseff <- ggplot(loss_df) +
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

# Weights -----------------------------------------------------------------

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

wt_loss <- weights_fun(model_model = loss_corr, 
                       model_raw = loss_raw)


corr_w <- as.data.frame(loss_corr$quantiles) %>%
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


# all ---------------------------------------------------------------------

modeled_losseff
ggsave(here('examples',
            'pictures',
            'original_R',
            'Fig4_modeledeffects_plot.pdf'),
       height =2,
       width = 3,
       units = "in")

weight_plot

ggsave(here('examples',
            'pictures',
            'original_R',
            'Fig4_weights_plot.pdf'),
       height =3,
       width = 3.75,
       units = "in")
