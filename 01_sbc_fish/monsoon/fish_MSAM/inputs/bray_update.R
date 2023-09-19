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

mod <- readRDS( file ="/scratch/atm234/sbc_fish/outputs/fish_MSAM_model.RDS")



# Update model to track bray ----------------------------------------------

parms <- c("bray")

mod2 <- update(mod,
               parameters.to.save = parms,
               n.iter = 100)

# Get samples from model --------------------------------------------------

#get the sims.list 
samples <- mod2$sims.list


# Select a site from the DF -----------------------------------------------

#site is the second dimension (columns)
#ABUR_1 is SiteID 1, so I'll just use it for now
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
  mutate(year = case_when(yearID == "V1" ~ 2000,
                          yearID == "V2" ~ 2001,
                          yearID == "V3" ~ 2002,
                          yearID == "V4" ~ 2003,
                          yearID == "V5" ~ 2004,
                          yearID == "V6" ~ 2005,
                          yearID == "V7" ~ 2006,
                          yearID == "V8" ~ 2007,
                          yearID == "V9" ~ 2008,
                          yearID == "V10" ~ 2009,
                          yearID == "V11" ~ 2010,
                          yearID == "V12" ~ 2011,
                          yearID == "V13" ~ 2012,
                          yearID == "V14" ~ 2013,
                          yearID == "V15" ~ 2014,
                          yearID == "V16" ~ 2015,
                          yearID == "V17" ~ 2016,
                          yearID == "V18" ~ 2017,
                          yearID == "V19" ~ 2018,
                          yearID == "V20" ~ 2019,
                          yearID == "V21" ~ 2020,
                          yearID == "V22" ~ 2021,
                          yearID == "V23" ~ 2022,
                          TRUE ~ NA_real_))

saveRDS(subsamples, file = "/scratch/atm234/sbc_fish/outputs/ABUR1_bray_samples.RDS")

mean <- mod2$mean$bray[1,]

sd <- mod2$sd$bray[1,]

year <- c(2000:2022)

summary <- as.data.frame(cbind(mean_bray = mean,
                         sd_bray = sd, 
                         year = year))

saveRDS(summary, file = "/scratch/atm234/sbc_fish/outputs/ABUR1_bray_summary.RDS")
