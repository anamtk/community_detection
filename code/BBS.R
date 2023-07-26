
# Load packages -----------------------------------------------------------


package.list <- c("here", "tidyverse")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load data ---------------------------------------------------------------

colo <- read.csv(here('data_raw',
                      "BBS",
                      "Colorad.csv"))

yrs <- colo %>%
  distinct(Route, Year) %>%
  group_by(Route) %>%
  tally()

colo2 <- colo %>%
  left_join(yrs, by = "Route") %>%
  filter(n > 19) %>%
  filter(!Year %in% c(2021, 2022))

yrID <- sort(unique(colo2$Year))

t1 <- colo2 %>%
  group_by(Route) %>%
  mutate(first_year = min(Year),
         last_year = max(Year)) %>%
  rowwise() %>%
  mutate(tot_years = last_year - first_year) %>%
  filter(n == tot_years) %>%
  ungroup()

t2 <- colo2 %>%
  distinct(Route, Year)
#only 7 have continuous data collection
t1 %>%
  distinct(Route) %>%
  tally()

