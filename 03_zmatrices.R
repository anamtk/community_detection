# Export z-matrices
# January 27, 2023
# Ana Miller-ter Kuile

#this script re-attaches metadata on species and years
# to the posterior z-matrices and exports them for downstream
# analyses

# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 

package.list <- c("here", "tidyverse")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

## And loading them

for(i in package.list){library(i, character.only = T)}
