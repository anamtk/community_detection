# preliminary community metrics stuff
# An Bui
# April 2023

# source ------------------------------------------------------------------

package.list <- c("here", "tidyverse", "rfishbase")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

## And loading them

for(i in package.list){library(i, character.only = T)}

source(here("code", "03_analyses", "community_metrics", "01_generating-list-cols.R"))

