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

#this is 100 random posterior samples of the "bray" object from the 
#model (what is the model predicting bray to be across specific
#iterations of the MCMC chains)
samples_sbc <- readRDS(here("05_visualizations",
                        "viz_data",
                        "ABUR1_bray_samples.RDS"))

#this is the posterior mean and standard deviation of bray
#from the model output
posterior_sbc <- readRDS(here("05_visualizations",
                          "viz_data",
                          "ABUR1_bray_summary.RDS"))

# ⊣ birds -----------------------------------------------------------------

#I output data for ABUR, transect 1 from the monthly fish survey
#data from 2002 - 2022

raw_knz <- readRDS(here("05_visualizations",
                         "viz_data",
                         "knz_004B_raw_bray.RDS"))

#this is 100 random posterior samples of the "bray" object from the 
#model (what is the model predicting bray to be across specific
#iterations of the MCMC chains)
samples_knz <- readRDS(here("05_visualizations",
                             "viz_data",
                             "knz_004B_bray_samples.RDS"))

#this is the posterior mean and standard deviation of bray
#from the model output
posterior_knz <- readRDS(here("05_visualizations",
                               "viz_data",
                               "knz_004B_bray_summary.RDS"))


# Grasshoppers ------------------------------------------------------------

raw_sev <- readRDS(here('05_visualizations',
                        'viz_data',
                        'sev_BOER_1_108_raw_bray.RDS'))

samples_sev <- readRDS(here('05_visualizations',
                            'viz_data',
                            'sev_BOER_1_108_bray_samples.RDS'))

posterior_sev <- readRDS(here('05_visualizations',
                              'viz_data',
                              'sev_BOER_1_108_bray_summary.RDS'))


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
  } else {
    warning("Check your arguments! You may have specified the wrong dataset.")
    return(NA)
  }
  
  if(dataset == "birds"){
    title = "KNZ birds"
  } else if(dataset == "fish") {
    title = "SBC fish"
  } else {
    warning("Check your arguments! You may have specified the wrong dataset.")
    return(NA)
  }
  
  if(dataset == "birds"){
    breaks = seq(from = 1982, to = 2010, by = 5)
  } else if(dataset == "fish") {
    breaks = seq(from = 2003, to = 2023, by = 4)
  } else {
    warning("Check your arguments! You may have specified the wrong dataset.")
    return(NA)
  }
  
  ggplot() +
    geom_line(data = samples_df,
              aes(x = year, y = bray, group = iter),
              color = modeled_col, alpha = 0.05) +
    theme(panel.grid = element_blank()) +
    geom_ribbon(data = posterior_df,
                aes(x = year, y = mean_bray, ymax = mean_bray + sd_bray, ymin = mean_bray - sd_bray), 
                fill = modeled_col, alpha = 0.3) +
    geom_line(data = posterior_df,
              aes(x = year, y = mean_bray),
              color = modeled_col, linewidth = 1) +
    geom_line(data = raw_df, 
              aes(x = year, y = raw_bray),
              color = observed_col, linewidth = 1) +
    labs(x = "Year", y = "Bray-Curtis dissimilarity",
         title = title) +
    scale_x_continuous(breaks = breaks) +
    theme(plot.title.position = "panel",
          plot.title = element_text(hjust = 0.5)) 
  
}

timeseries_sbc <- timeseries_function(dataset = "fish") +
  coord_cartesian(xlim = c(2003, 2022))
timeseries_knz <- timeseries_function(dataset = "birds") +
  coord_cartesian(xlim = c(1982, 2009))

timeseries_together <- timeseries_sbc | timeseries_knz
timeseries_together


