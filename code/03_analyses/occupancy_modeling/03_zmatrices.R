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


# Load z-matrices ---------------------------------------------------------

z <- readRDS(here('data_outputs',
                  'monsoon_outputs',
                  'fish_MSOM_1_26_matrices.RDS'))


# Load raw matrix with info on rows and column names ----------------------

raw <- readRDS(here("data_outputs",
                     "community_matrices",
                     "fish_AQUE1_raw_matrix.RDS"))

species <- rownames(raw)
years <- colnames(raw)

# Add row and column names to z-matrices ----------------------------------

rownames(z[[1]]) <- species
rownames(z[[2]]) <- species
rownames(z[[3]]) <- species
colnames(z[[1]]) <- years
colnames(z[[2]]) <- years
colnames(z[[3]]) <- years


# Export z matrices with metadata -----------------------------------------

saveRDS(z, here("data_outputs",
                "community_matrices",
                "fish_AQUE1_quantile_matrices.RDS"))
