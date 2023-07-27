# look at posteriors for MSOM
# Ana Miller-ter Kuile
# January 23, 2023

# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 

package.list <- c("here", "tidyverse")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

## And loading them

for(i in package.list){library(i, character.only = T)}

# Load posterior data -----------------------------------------------------

summary <- readRDS(here("monsoon",
                             "outputs",
                             "sum_detection_covs.RDS"))

# Look at distributions of posteriors -------------------------------------

View(summary$quantiles)

