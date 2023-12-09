#Loading and description of raw vs. corrected graph dat objects
#Ana Miller-ter Kuile
#September 20, 2023

#this script loads and describes the three datasets generated
#to create the raw vs. corrected Bray-Curtis figure 

# Load packages -----------------------------------------------------------


package.list <- c("here", "tidyverse")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load data ---------------------------------------------------------------

# ⊣ sbc ------------------------------------------------------------------

#I output data for ABUR, transect 1 from the monthly sbc survey
#data from 2002 - 2022

raw_sbc <- readRDS(here("05_visualizations",
                    "viz_data",
                    "sbc_ABUR1_raw_bray.RDS"))

#rename to merge
raw_sbc <- raw_sbc %>%
  rename(raw_diss = raw_bray_all) 

#this is 100 random posterior samples of the "bray" object from the 
#model (what is the model predicting bray to be across specific
#iterations of the MCMC chains)
samples_sbc <- readRDS(here("05_visualizations",
                        "viz_data",
                        "ABUR1_bray_samples.RDS"))

samples_sbc <- samples_sbc %>%
  rename(diss = bray)


#this is the posterior mean and standard deviation of bray
#from the model output
posterior_sbc <- readRDS(here("05_visualizations",
                          "viz_data",
                          "ABUR1_bray_summary.RDS"))

posterior_sbc <- posterior_sbc %>%
  rename(mean = mean_bray,
         sd = sd_bray)

# ⊣ birds -----------------------------------------------------------------

raw_knz <- readRDS(here("05_visualizations",
                         "viz_data",
                         "knz_N04D_raw_bray.RDS"))

raw_knz <- raw_knz %>%
  rename(raw_diss = raw_bray_all)

#this is 100 random posterior samples of the "bray" object from the 
#model (what is the model predicting bray to be across specific
#iterations of the MCMC chains)
samples_knz <- readRDS(here("05_visualizations",
                             "viz_data",
                             "knz_N04D_bray_samples.RDS"))

samples_knz <- samples_knz %>%
  rename(diss = bray)

#this is the posterior mean and standard deviation of bray
#from the model output
posterior_knz <- readRDS(here("05_visualizations",
                               "viz_data",
                               "knz_N04D_bray_summary.RDS"))

posterior_knz <- posterior_knz %>%
  rename(mean = mean_bray,
         sd = sd_bray)

# ⊣ grasshoppers ----------------------------------------------------------

raw_sev <- readRDS(here('05_visualizations',
                        'viz_data',
                        'sev_BOER_1_108_raw_bray.RDS'))

raw_sev <- raw_sev %>%
  rename(raw_diss = raw_bray_all)

samples_sev <- readRDS(here('05_visualizations',
                            'viz_data',
                            'sev_BOER_1_108_bray_samples.RDS'))

samples_sev <- samples_sev %>%
  rename(diss = bray)

posterior_sev <- readRDS(here('05_visualizations',
                              'viz_data',
                              'sev_BOER_1_108_bray_summary.RDS'))

posterior_sev <- posterior_sev %>%
  rename(mean = mean_bray,
         sd = sd_bray)


# Plants ------------------------------------------------------------------

nps_ids <- read.csv(here('04_nps_plants',
                         'data_outputs',
                         'metadata',
                         'site_year_IDs.csv'))

raw_nps <- readRDS(here('05_visualizations',
                        'viz_data',
                        'nps_S02_B_3_raw_jaccard.RDS'))

raw_nps <- raw_nps %>%
  rename(raw_diss = turnover_all)


samples_nps <- readRDS(here('05_visualizations',
                            'viz_data',
                            'nps_Jaccard_samples_S02_B_3.RDS'))

years <- nps_ids %>%
  unite(c("Plot", "Transect", "Quadrat"),
        col = "plot_trans_quad",
        sep = "_",
        remove = F) %>%
  distinct(plot_trans_quad, EventYear) %>%
  group_by(plot_trans_quad) %>%
  arrange(EventYear) %>%
  mutate(yrNum = 1:n()) %>%
  filter(plot_trans_quad == "S02_B_3") %>%
  ungroup() %>%
  dplyr::select(-plot_trans_quad)

samples_nps <- samples_nps %>%
  left_join(years, by = "yrNum") %>%
  rename('iter' = 'groups',
         'yearID' = 'yrNum',
         'diss' = 'turnover',
         'year' = 'EventYear') %>%
  dplyr::select(-plot_trans_quad)

posterior_nps <- readRDS(here('05_visualizations',
                              'viz_data',
                              'nps_Jaccard_summary_S02_B_3.RDS'))

posterior_nps <- posterior_nps %>%
  rename(year = EventYear) %>%
  dplyr::select(-plot_trans_quad)

# plot --------------------------------------------------------------------

modeled_col <- "#E88C23"
observed_col <- "#438AA8"

timeseries_function <- function(dataset) {
  
  if(dataset == "birds"){
    samples_df <- samples_knz
    posterior_df <- posterior_knz
    raw_df <- raw_knz
  } else if(dataset == "fish") {
    samples_df <- samples_sbc
    posterior_df <- posterior_sbc
    raw_df <- raw_sbc
  } else if(dataset == "grasshoppers") {
    samples_df <- samples_sev
    posterior_df <- posterior_sev
    raw_df <- raw_sev
  } else if(dataset == "plants") {
    samples_df <- samples_nps
    posterior_df <- posterior_nps
    raw_df <- raw_nps
  } else {
    warning("Check your arguments! You may have specified the wrong dataset.")
    return(NA)
  }
  
  if(dataset == "birds"){
    title = "(b) KNZ birds"
  } else if(dataset == "fish") {
    title = "(a) SBC fish"
  } else if(dataset == "grasshoppers") {
    title = "(c) SEV grasshoppers"
  } else if(dataset == "plants") {
    title = "(d) PFNP plants"
  } else {
    warning("Check your arguments! You may have specified the wrong dataset.")
    return(NA)
  }
  
  if(dataset == "birds"){
    breaks = seq(from = 1982, to = 2010, by = 5)
  } else if(dataset == "fish") {
    breaks = seq(from = 2003, to = 2023, by = 4)
  } else if(dataset == "grasshoppers") {
    breaks = seq(from = 1993, to = 2019, by = 5)
  } else if(dataset == "plants") {
    breaks = seq(from = 2008, to = 2024, by = 3)
  } else {
    warning("Check your arguments! You may have specified the wrong dataset.")
    return(NA)
  }
  
  if(dataset == "birds"){
    xlab = ""
  } else if(dataset == "fish") {
    xlab = ""
  } else if(dataset == "grasshoppers") {
    xlab = ""
  } else if(dataset == "plants") {
    xlab = "Year"
  } else {
    warning("Check your arguments! You may have specified the wrong dataset.")
    return(NA)
  }
  
  ggplot() +
    geom_line(data = samples_df,
              aes(x = year, y = diss, group = iter),
              color = modeled_col, alpha = 0.05) +
    theme(panel.grid = element_blank()) +
    geom_ribbon(data = posterior_df,
                aes(x = year, y = mean, ymax = mean + sd, ymin = mean - sd), 
                fill = modeled_col, alpha = 0.3) +
    geom_line(data = posterior_df,
              aes(x = year, y = mean),
              color = modeled_col, linewidth = 1) +
    geom_line(data = raw_df, 
              aes(x = year, y = raw_diss),
              color = observed_col, linewidth = 1) +
    labs(x = xlab, y = "Dissimilarity",
         title = title) +
    scale_x_continuous(breaks = breaks) +
    theme(legend.position = "none",
          plot.title.position = "plot",
          panel.grid = element_blank(),
          axis.text = element_text(size = 6),
          axis.title = element_text(size = 7),
          plot.title = element_text(size = 8)) 
  
}

timeseries_sbc <- timeseries_function(dataset = "fish") +
  coord_cartesian(xlim = c(2003, 2022)) 
timeseries_knz <- timeseries_function(dataset = "birds") +
  coord_cartesian(xlim = c(1982, 2009))
timeseries_sev <- timeseries_function(dataset = "grasshoppers") +
  coord_cartesian(xlim = c(1993, 2019))
timeseries_nps <- timeseries_function(dataset = "plants") +
  coord_cartesian(xlim = c(2008, 2022))

timeseries_together <- timeseries_sbc / timeseries_knz / timeseries_sev / timeseries_nps
timeseries_together 

ggsave(plot = timeseries_together,
       filename = here("pictures",
                       "detection_models",
                       "observed_modeled_timeseries.jpg"),
       height = 14,
       width = 12,
       units = "cm")
