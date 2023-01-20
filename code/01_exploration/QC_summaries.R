# Ana Miller-ter Kuile
# January 13, 2023

# this script explores the algae cover dataset, examining the species
# that are present, and looking at their distribution of covers
# and maybe generating output data for the model and/or average
# covers per species

# Load packages -----------------------------------------------------------


package.list <- c("here", "tidyverse", "patchwork")

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
#Stephanocystis osmundacea - CH and CYOS
#Macrocystis pyrifera - DMH and MH
#Amphiroa beauvoisii - AMZO, UEC
benthic1 <- benthic %>%
  filter(!str_detect(SCIENTIFIC_NAME, "Unidentified|Unidentifiable")) %>%
  mutate(SP_CODE = case_when(SCIENTIFIC_NAME == "Stephanocystis osmundacea" ~ "CH",
                             SCIENTIFIC_NAME == "Macrocystis pyrifera" ~ "DMH",
                             SCIENTIFIC_NAME == "Amphiroa beauvoisii" ~ 'AMZO',
                             TRUE ~ SP_CODE))


benthic1 %>%
  filter(PERCENT_COVER > 0) %>%
  group_by(YEAR, MONTH, SITE, TRANSECT, QUAD, SIDE) %>%
  tally() %>%
  ggplot(aes(x = n)) + 
  geom_histogram()

benthic1 %>%
  filter(PERCENT_COVER > 0) %>%
  group_by(SCIENTIFIC_NAME) %>%
  summarise(mean = mean(PERCENT_COVER),
            sd = sd(PERCENT_COVER),
            total = n(),
            se = sd/sqrt(total)) %>%
  ggplot(aes(x = reorder(SCIENTIFIC_NAME, mean), y = mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2) +
  coord_flip()

#things labeled "spp" - are these one species? 
# what about "unidentified things?
# how to get to species level, or how should we split things up - based 
## on what you ahve done before?
# things with multiple scientific names? are these synonyms?
# is there a dataset of the average individual body size for these species?
## just trying to brainstorm how to use the percent cover type data as a 
## covariate for detection probability
