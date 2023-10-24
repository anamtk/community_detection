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

#I output data for ABUR, transect 1 from the monthly fish survey
#data from 2002 - 2022

raw <- readRDS(here("05_visualizations",
                    "viz_data",
                    "sbc_ABUR1_raw_bray.RDS"))

#this is 100 random posterior samples of the "bray" object from the 
#model (what is the model predicting bray to be across specific
#iterations of the MCMC chains)
samples <- readRDS(here("05_visualizations",
                        "viz_data",
                        "ABUR1_bray_samples.RDS"))

#this is the posterior mean and standard deviation of bray
#from the model output
posterior <- readRDS(here("05_visualizations",
                          "viz_data",
                          "ABUR1_bray_summary.RDS"))


# plot --------------------------------------------------------------------

modeled_col <- "#E88C23"
observed_col <- "#438AA8"

ggplot() +
  geom_line(data = samples,
            aes(x = year, y = bray, group = iter),
            color = modeled_col, alpha = 0.05) +
  theme(panel.grid = element_blank()) +
  geom_ribbon(data = posterior,
              aes(x = year, y = mean_bray, ymax = mean_bray + sd_bray, ymin = mean_bray - sd_bray), 
              fill = modeled_col, alpha = 0.3) +
  geom_line(data = posterior,
            aes(x = year, y = mean_bray),
            color = modeled_col, linewidth = 1) +
  geom_line(data = raw, 
            aes(x = year, y = raw_bray),
            color = observed_col, linewidth = 1)


# Figure thoughts ---------------------------------------------------------

#I really like your original figure. I wonder if ew could add in 
#what you already had (raw vs 100 iterations) and then the 
#mean and standard deviation of the posterior as another line 
#(with error bars or a shaded ribbon?) 

#anyway, those are my thoughts - I'm sure you'll create something
#beautiful!

#Ideally, we will be generating some version of this figure for 
#all the other datasets as well, so if you can create the code to
#be somewhat transferrable among datasets with slightly different
#structures, that could be great!

