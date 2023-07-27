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

summaries <- readRDS(here("monsoon",
                          "outputs",
                          "fish_community_summaries.RDS"))

ids <- read.csv(here("data_outputs",
                      "metadata",
                      "site_year_IDs.csv"))


# Pull out medians --------------------------------------------------------

meds <- as.data.frame(summaries$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(parm != "deviance") %>%
  separate(parm,
           into = c("parm", "site_year"),
           sep = "\\[") %>%
  separate(site_year,
           into = c("siteID", "yearID"),
           sep = ",") %>%
  mutate(yearID = str_sub(yearID, 1, end = -2L)) %>%
  rename('median' = '50%') %>%
  dplyr::select(parm, siteID, yearID, median) %>%
  pivot_wider(names_from = parm,
              values_from = median) %>%
  mutate(siteID = as.numeric(siteID),
         yearID = as.numeric(yearID)) %>%
  left_join(ids, by = c("siteID", "yearID" = "yrID"))

write.csv(meds, here("data_outputs",
                     "community_stability",
                     "corrected_stability_metrics.csv"))
