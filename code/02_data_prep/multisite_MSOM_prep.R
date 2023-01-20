# Ana Miller-ter Kuile
# January 17, 2023



# Load packages -----------------------------------------------------------


package.list <- c("here", "tidyverse")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load data ---------------------------------------------------------------

benthic <- read.csv(here("data_raw", 
                         "Annual_Cover_All_Years_20220809.csv"))



# Explore benthic dataset -------------------------------------------------

colnames(benthic)

#2 species codes for
#Stephanocystis osmundacea - CH and CYOS - this is holdfast vs. young
#Macrocystis pyrifera - DMH and MH - this is holdfast vs. young
#Amphiroa beauvoisii - AMZO, UEC - this is two spp code with no explanation...

#remove unidentified species from the group
benthic1 <- benthic %>%
  filter(!str_detect(SCIENTIFIC_NAME, "Unidentified|Unidentifiable")) %>%
  #one species code is confusingly "NA" so have to tell computer
  # that isn't missing data
  mutate(SP_CODE = case_when(is.na(SP_CODE) ~ "NA",
                             TRUE ~ SP_CODE)) %>%
  group_by(SCIENTIFIC_NAME, SITE, TRANSECT, QUAD, SIDE, YEAR) %>%
  #some species have more than one code, so combining here (might
  # change this appraoch later)
  summarise(PERCENT_COVER1 = sum(PERCENT_COVER)) %>%
  ungroup()

benthic2 <- benthic1 %>%
  #get occupancy from percent cover
  mutate(OCC = case_when(PERCENT_COVER1 > 0 ~ 1,
                         PERCENT_COVER1 == 0 ~ 0,
                         TRUE ~ NA_real_)) %>%
  #get unique side_quad combos
  unite(SIDE_QUAD, c("QUAD", "SIDE"), sep = "_") %>%
  #get unique site_transect combos
  unite(SITE_TRANSECT, c("SITE", "TRANSECT"), sep = "_") %>%
  #select variables of interest
  dplyr::select(YEAR,SIDE_QUAD, SITE_TRANSECT, SCIENTIFIC_NAME, OCC) %>%
  #remove missing occupancy values now 
  #might want to set to zero? Think about this later down the road
  filter(!is.na(OCC)) 

#this gets a year number - maybe not necessary -we'll see
benthic_meta <- benthic2 %>%
  distinct(SITE_TRANSECT, YEAR) %>%
  group_by(SITE_TRANSECT) %>%
  arrange(SITE_TRANSECT, YEAR) %>%
  mutate(YEAR_NUM = 1:n()) %>%
  ungroup()


t <- benthic2 %>%
  left_join(benthic_meta, by = c("SITE_TRANSECT", "YEAR")) %>%
  filter((YEAR_NUM == 1 & SIDE_QUAD == '20_I')) %>%
  dplyr::select(-YEAR_NUM, -SIDE_QUAD, -YEAR) %>%
  pivot_wider(names_from = SITE_TRANSECT,
              values_from = OCC) %>% 
  column_to_rownames(var = "SCIENTIFIC_NAME")

t <- benthic2 %>%
  left_join(benthic_meta, by = c("SITE_TRANSECT", "YEAR")) %>%
  filter((YEAR_NUM == 13 & SIDE_QUAD == '20_I')) %>%
  dplyr::select(-YEAR_NUM, -SIDE_QUAD, -YEAR) %>%
  pivot_wider(names_from = SITE_TRANSECT,
              values_from = OCC) %>% 
  column_to_rownames(var = "SCIENTIFIC_NAME")

