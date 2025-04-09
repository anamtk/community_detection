#model values results
#Ana Miller-ter Kuile
#February 20, 2025

#summariseing results from modeled/observed values of all datasets

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

grasshopper<- read.csv(here('examples',
                     'data_output',
                     'grasshopper_stability',
                     'tidy_data',
                     "grasshopper_betareg_tidydata.csv"))


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

plant <- read.csv(here('examples',
                          'data_output',
                          'plant_turnover',
                          'tidy_data',
                          'plant_betareg_tidydata.csv'))


# Manipulate and plot -----------------------------------------------------

#how do raw/correcetd compare?
bird_rao <- bird_all %>%
  dplyr::select(TransID, yrID, Q_mean, Q) %>%
  rename(modeled = Q_mean,
         observed = Q) %>%
  pivot_longer(modeled:observed,
               names_to = "type",
               values_to = "Q")

(bird_plot <- ggplot(bird_rao, aes(x = type,
                                   fill = type,
                                   y = Q)) +
    geom_boxplot()+
    scale_fill_manual(values = c(modeled_col, observed_col),
                      labels = c("Modeled", "Empirical")) +
    theme(legend.position = "none") +
    labs(x = "Data type",
         y = "Rao's quadratic entropy"))

grasshopper2 <- grasshopper %>%
  dplyr::select(mean_bray, observed_all) %>%
  rename('modeled' = mean_bray,
         observed = observed_all) %>%
  pivot_longer(modeled:observed,
               names_to = "type",
               values_to = 'bray')

(hop_plot <- ggplot(grasshopper2, aes(x = type,
                                   fill = type,
                                   y = bray)) +
    geom_boxplot()+
    scale_fill_manual(values = c(modeled_col, observed_col),
                      labels = c("Modeled", "Empirical")) +
    theme(legend.position = "none") +
    labs(x = "Data type",
         y = "Bray-Curtis temporal dissimilarity"))

plant_loss <- plant %>%
  dplyr::select(loss, mean_loss) %>%
  rename(modeled = mean_loss,
         observed = loss) %>%
  pivot_longer(modeled:observed,
               names_to = "type",
               values_to = "loss")

(plant_plot_box <- ggplot(plant_loss, aes(x = type, y = loss, fill = type)) +
    geom_boxplot()+
    scale_fill_manual(values = c(modeled_col, observed_col),
                      labels = c("Modeled", "Empirical"))+
    labs(x = "Data type",
         y = "Species losses") )


# all plot ----------------------------------------------------------------

bird_plot +plant_plot_box+hop_plot +
  plot_annotation(tag_levels = "a",
                  tag_suffix = ")") +
  plot_layout(guides = "collect")

ggsave(here('examples',
            'pictures',
            'original_R',
            'Fig2b_diversity_metrics.pdf'),
       height = 3, 
       width = 10,
       units = "in")

