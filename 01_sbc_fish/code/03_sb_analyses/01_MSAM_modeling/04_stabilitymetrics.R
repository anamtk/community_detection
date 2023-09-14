#Stability by community-year
#Ana Miller-ter Kuile
#May 15, 2023

#this script pulls out the median stability metrics per community/year
#adn gives them their original data IDs to run follow-up analyses

# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 

package.list <- c("here", "tidyverse")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

## And loading them

for(i in package.list){library(i, character.only = T)}


# Load data ---------------------------------------------------------------

summaries <- readRDS(here('sbc_fish',
                          "monsoon",
                          "outputs",
                          "fish_community_summaries.RDS"))

sum <- summary(mod$samples)

summaries <- as.data.frame(sum$statistics) %>%
  rownames_to_column(var = "parameter") %>%
  dplyr::select(parameter, Mean, SD) %>%
  filter(parameter != 'deviance') %>%
  separate(parameter,
           into = c("parm", "type"),
           sep = '\\[') %>%
  mutate(type = str_sub(type, start =1, end = nchar(type)-1)) %>%
  separate(type,
           into = c("siteID", "yrID"),
           sep = ',') %>%
  mutate(siteID = as.numeric(siteID),
         yrID = as.numeric(yrID))

ids <- read.csv(here('sbc_fish',
                     "data_outputs",
                      "metadata",
                      "site_year_IDs.csv"))


# Pull out medians --------------------------------------------------------

dat2 <- summaries %>%
  left_join(ids, by = c("siteID", "yrID"))

write.csv(dat2, here('sbc_fish',
                     "data_outputs",
                     'SAM',
                     "data_prep",
                     "corrected_stability_metrics.csv"))
