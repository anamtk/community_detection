#Get Bray-Curtis out of converged model
#Ana Miller-ter Kuile
#September 19, 2023

#this script pulls out and sumamrises Bray Curtis for visualizations
#purposes for An's project


# Load packages ---------------------------------------------------------------
Sys.time()


# Load packages,
package.list <- c('jagsUI',"coda",'dplyr', 
                   'stringr','magrittr',
                   'tidyr','ggplot2',
                   'tibble') 


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load model --------------------------------------------------------------

mod <- readRDS( file ="/scratch/atm234/sev_hoppers/outputs/sevMSAMnocov_model2.RDS")



# Update model to track bray ----------------------------------------------

parms <- c("bray")

mod2 <- update(mod,
               parameters.to.save = parms,
               n.iter = 4000)


# I need mean and SD for all sites ----------------------------------------

sum <- summary(mod2$samples)

#matrix of site x year:
stats <- sum$statistics

saveRDS(stats, 
        file = "/scratch/atm234/sev_hoppers/outputs/sev_bray_meanSD.RDS")

# Get samples from model --------------------------------------------------

#get the sims.list 
samples <- mod2$sims.list


# Select a site from the DF -----------------------------------------------

#site is the second dimension (columns)
# BOER_1_108 is SiteID 1, so I'll just use it for now
site1 <- as.data.frame(samples$bray[,1,])

#how many samples to select from the DF
nsamps <- 100

#subset a random set of nsamps samples from either samples or sims.list
subsamples <- slice_sample(.data = site1, 
                           n = nsamps,
                           replace= FALSE) %>%
  #give each sample an ID
  mutate(iter = 1:n()) %>%
  #pivot longer
  pivot_longer(V1:V23,
               names_to = "yearID",
               values_to = "bray") %>%
  #give years sensical IDs
  mutate(year = case_when(yearID == "V1" ~ 1992,
                          yearID == "V2" ~ 1993,
                          yearID == "V3" ~ 1994,
                          yearID == "V4" ~ 1995,
                          yearID == "V5" ~ 1996,
                          yearID == "V6" ~ 1997,
                          yearID == "V7" ~ 1998,
                          yearID == "V8" ~ 1999,
                          yearID == "V9" ~ 2000,
                          yearID == "V10" ~ 2001,
                          yearID == "V11" ~ 2002,
                          yearID == "V12" ~ 2003,
                          yearID == "V13" ~ 2004,
                          yearID == "V14" ~ 2005,
                          yearID == "V15" ~ 2006,
                          yearID == "V16" ~ 2007,
                          yearID == "V17" ~ 2008,
                          yearID == "V18" ~ 2009,
                          yearID == "V19" ~ 2010,
                          yearID == "V20" ~ 2011,
                          yearID == "V21" ~ 2012,
                          yearID == "V22" ~ 2013,
                          yearID == "V23" ~ 2014,
                          yearID == "V23" ~ 2015,
                          yearID == "V23" ~ 2016,
                          yearID == "V23" ~ 2017,
                          yearID == "V23" ~ 2018,
                          yearID == "V23" ~ 2019,
                          TRUE ~ NA_real_))

saveRDS(subsamples, file = "/scratch/atm234/sev_hoppers/outputs/sev_BOER_1_108_bray_samples.RDS")

mean <- mod2$mean$bray[1,]

sd <- mod2$sd$bray[1,]

year <- c(1992:2019)

summary <- as.data.frame(cbind(mean_bray = mean,
                         sd_bray = sd, 
                         year = year))

saveRDS(summary, file = "/scratch/atm234/sev_hoppers/outputs/sev_BOER_1_108_bray_summary.RDS")
