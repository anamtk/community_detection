# Export z-matrices
# January 27, 2023
# Ana Miller-ter Kuile

#this script re-attaches metadata on species and years
# to the posterior z-matrices and exports them for downstream
# analyses

# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 

package.list <- c("here", "tidyverse",
                  "coda")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

## And loading them

for(i in package.list){library(i, character.only = T)}


# Load z-matrices ---------------------------------------------------------

z <- readRDS(here('sbc_fish',
                  "data_outputs",
                  "community_matrices",
                  "fish_MSOM_community_matrices.RDS"))

#these are in a WILD format where rows = iteration number,
# columns are of the form 'z[x,y]', where the 1-0 value
# corresponds to the row and column value for that matrix in that
# iteration. x = species, y = year

# Combine all Z's ---------------------------------------------------------

#pull all three 500 iteration chains out of the weird mcmc object
a <- as.matrix(z[[1]])
b <- as.matrix(z[[2]])
c <- as.matrix(z[[3]])

#make this a dataframe bound by rows
z_all <- as.data.frame(rbind(a,b,c))

#subset just a few for now because this is a LOT - can change the n 
# later if we want - I would probably aim for ~1000+ for any final
# analyses
z_all2 <- slice_sample(z_all, n = 100)


# Get all z matrix into long format ---------------------------------------

z_all3 <- z_all2 %>%
  #remove the deviance column
  dplyr::select(-deviance) %>%
  #make the rownames the "iteration"
  rownames_to_column(var = "iteration") %>%
  #pivot longer everything but the iteration value
  pivot_longer(2:last_col(),
               #each new line is now a species-year combo for that iteration
               names_to = 'species_year',
               #and occupancy gets its own column
               values_to = "occupancy") %>%
  #get all the weird characters out of each row cell value for this
  #species-year column (remove the Z and the [])
  mutate(species_year = str_sub(species_year, 
                                start = 3, 
                                end = str_length(species_year)-1)) %>%
  #separate that column into a species ID and a year ID
  separate(species_year, into = c("species_num", "year_num"), sep = ",") %>%
  #make these numeric so we can bind them together later
  mutate(species_num = as.numeric(species_num),
         year_num = as.numeric(year_num))

# Load raw matrix with info on rows and column names ----------------------

raw <- readRDS(here("data_outputs",
                     "community_matrices",
                     "fish_AQUE1_raw_matrix.RDS"))

#get species names linked to specifci numerical id that corresponds to above
#matrix
species <- as.data.frame(rownames(raw)) %>%
  mutate(species_num = 1:n()) %>%
  rename("species" = "rownames(raw)")

#get year ids linked to specific numerical id that corresponds to above matrix
years <- as.data.frame(colnames(raw)) %>%
  mutate(year_num = 1:n()) %>%
  rename("year" = "colnames(raw)")

# Add year and species values to the z-matrices ---------------------------

#add the year and species IDs to the long dataframe
z_all4 <- z_all3 %>%
  left_join(species, by = "species_num") %>%
  left_join(years, by = "year_num")


# Export z matrices with metadata -----------------------------------------

write.csv(z_all4, here("data_outputs",
                "community_matrices",
                "fish_AQUE1_z_matrices.csv"))
