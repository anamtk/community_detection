# preliminary trait stuff
# An Bui
# July 2023

# source ------------------------------------------------------------------

package.list <- c("here", "tidyverse", "rfishbase")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

## And loading them

for(i in package.list){library(i, character.only = T)}

# pulling traits from FishBase --------------------------------------------

# creating list of genus and species for fish
sbc_spp <- read_rds(here("data_outputs", "metadata", "species_metadata.RDS")) %>%
  # filtering out scientific names where no fish were detected and unidentified species
  filter(!(scientific_name %in% c(".",
                                  "Unidentified Flatfish spp.", 
                                  "Unidentified Surfperch spp.",
                                  # kelp fish
                                  "Unidentified YOY Sebastes",
                                  # sculpins
                                  "Unidentified Cottidae",
                                  # pipefish genus
                                  "Syngnathus spp.",
                                  # kelp fish genus, include all 3 species?
                                  "Gibbonsia"))) %>% 
  pull(scientific_name)

# check to make sure that names are in FishBase
validate_names(sbc_spp)

# pulling out traits for SBC species
# climate vulnerability
sbc_spp_traits <- fb_tbl("species") %>% 
  mutate(sci_name = paste(Genus, Species)) %>%
  filter(sci_name %in% sbc_spp)

# write_rds(x = sbc_spp_traits, file = here("data_outputs", "traits", "sbc_spp_traits.RDS"))

# trophic level
sbc_spp_diet <- diet(species_list = sbc_spp)

# offspring information
sbc_spp_fecundity <- fecundity(species_list = sbc_spp)

# age and length at maturity
sbc_spp_maturity <- maturity(species_list = sbc_spp)

# feeding type
sbc_spp_ecology <- ecology(species_list = sbc_spp)

# descriptions of swimming movements
sbc_spp_swimming <- swimming(species_list = sbc_spp)

