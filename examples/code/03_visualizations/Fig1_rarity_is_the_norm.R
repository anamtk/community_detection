
# Load packages -----------------------------------------------------------


package.list <- c("here", "tidyverse", 'ratdat',
                  'patchwork')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

theme_set(theme_bw())

# Load data ---------------------------------------------------------------

#get link to this dataset from sbc lter
fish <- read.csv(here('examples',
                      'data_raw',
                      'other_datasets',
                      'sbc_fish_data.csv'))

birds <- read.csv(here('examples',
                       'data_output',
                       'bird_fundiv',
                       'tidy_data',
                       'bird_msam_tidy_data.csv'))

hoppers <- read.csv(here('examples',
                         'data_output',
                         'grasshopper_stability',
                         'tidy_data',
                         'grasshopper_msam_tidy_data.csv'))

plants <- read.csv(here('examples',
                        'data_output',
                        'plant_turnover',
                        'tidy_data',
                        'plant_msom_tidy_data.csv'))

rodents <- ratdat::surveys

#from Phoenix LTER:
#https://globalfutures.asu.edu/caplter/
#https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-cap.627.9
herps <- read.csv(here('examples',
                       'data_raw',
                       'other_datasets',
                       '627_herp_survey_observations.csv'))

#ants
#https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.1682#support-information-section
ants <- read.csv(here('examples',
                      'data_raw',
                      'other_datasets',
                      'Gibb et al. Ecology - observations data.csv'))
# Summarise rarity --------------------------------------------------------

#two ways to think about rarity:
#low abundance overall 
#low frequency of occurrence across space/time

# Frequency through space time --------------------------------------------

freq_fun <- function(df, yearID, siteID, repID, countID, speciesID, title){
  
  total <- df %>%
    distinct({{yearID}}, {{siteID}}, {{repID}}) %>%
    summarise(total = n()) %>%
    ungroup() %>%
    dplyr::select(total) %>%
    as_vector()
  
  df2 <- df %>%
    filter({{countID}} > 0) %>%
    group_by({{speciesID}}) %>%
    tally() %>%
    ungroup() %>%
    rowwise() %>%
    summarise(freq= n/total) %>%
    ungroup() %>%
    mutate(dataset = title)# %>%
    #dplyr::select(-{{speciesID}})
  
  plot <- ggplot(df2) +
    geom_histogram(aes(x = n)) +
    labs(x = "Frequency of observation",
         y = "Number of species",
         title = title) 
  
  return(df2)
    
}

fish_freq <- freq_fun(df = fish,
                      yearID = YEAR,
                      siteID = SITE_TRANS, 
                      repID = REP,
                      countID = COUNT,
                      speciesID = SP_CODE,
                      title = "fish")

bird_freq <- freq_fun(df = birds,
                      yearID = yrID,
                      siteID = TransID,
                      repID = REP,
                      speciesID = SpecID,
                      countID = NOBS,
                      title = "birds")

hopper_freq <- freq_fun(df = hoppers,
                        yearID = yrID,
                        siteID = siteID,
                        repID = rep,
                        speciesID = speciesID,
                        countID = CNT,
                        title = "grasshoppers")

plant_freq <- freq_fun(df = plants,
                       yearID = yrID,
                       siteID = quadnum,
                       repID = REP,
                       speciesID = SpecID,
                       countID = presence,
                       title = "plants")

rodent_total <- rodents %>%
  distinct(month, day, year, plot_id) %>%
  tally() %>%
  ungroup() %>%
  dplyr::select(n) %>%
  as_vector()

rodent_freq <- rodents %>%
  group_by(month, day, year, plot_id, species_id) %>%
  summarise(count = n()) %>%
  ungroup()%>%
  group_by(species_id) %>%
  tally() %>%
  ungroup() %>%
  rowwise() %>%
  summarise(freq = n/rodent_total) %>%
  ungroup() %>%
  mutate(dataset = 'rodents')
  
herp_total <- herps %>%
  distinct(reach, transect, location, 
           observation_date) %>%
  tally() %>%
  ungroup() %>%
  dplyr::select(n) %>%
  as_vector()
  
herp_freq <- herps %>%
  group_by(reach, transect, location, 
           observation_date, common_name) %>%
  summarise(count = n()) %>%
  ungroup()%>%
  group_by(common_name) %>%
  tally() %>%
  ungroup() %>%
  rowwise() %>%
  summarise(freq = n/herp_total) %>%
  ungroup() %>%
  mutate(dataset = 'herps')
  

freq_df <- bind_rows(fish_freq, hopper_freq,
                     bird_freq, plant_freq,
                     herp_freq, rodent_freq)

b <- ggplot(freq_df) + 
  geom_histogram(aes(x = freq, fill = dataset),
                 position = position_dodge()) +
  scale_x_sqrt() +
  facet_grid(dataset~., scales = "free_y") +
  labs(x = "Frequency of observation \n(proportion of surveys)",
       y = "Number of species") +
  theme(legend.position = 'none',
        strip.background = element_rect(fill = "white"))

ggplot(freq_df) + 
  geom_density(aes(x = freq, fill = dataset), alpha = 0.3) +
  scale_x_sqrt() +
  labs(x = "Frequency of observation \n(proportion of surveys)",
       y = "Number of species") +
  theme(legend.position = 'none',
        strip.background = element_rect(fill = "white"))

# All surveyes ------------------------------------------------------------

fish_all <- fish %>%
  dplyr::select(COUNT) %>%
  mutate(dataset = "fish")

hop_all <- hoppers %>%
  dplyr::select(CNT) %>%
  rename(COUNT = CNT) %>%
  mutate(dataset = 'grasshoppers')

bird_all <- birds %>%
  dplyr::select(NOBS) %>%
  rename(COUNT = NOBS) %>%
  mutate(dataset = "birds")

rodent_all <- rodents %>%
  group_by(month, day, year, plot_id, species_id) %>%
  summarise(count = n()) %>%
  ungroup()%>%
  mutate(species_id = as.factor(species_id)) %>%
  group_by(month, day, year, plot_id) %>%
  complete(species_id) %>%
  replace_na(list(count = 0)) %>% 
  dplyr::select(count) %>%
  rename(COUNT = count) %>%
  mutate(dataset = "rodents")

herps_all <- herps %>%
  filter(!is.na(common_name)) %>%
  group_by(reach, transect, location, 
           observation_date,common_name) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  mutate(common_name = as.factor(common_name)) %>%
  group_by(reach, transect, location, 
           observation_date)%>%
  complete(common_name) %>%
  ungroup()  %>%
  replace_na(list(count = 0)) %>%
  dplyr::select(count) %>%
  rename(COUNT = count) %>%
  mutate(dataset = 'herps')


ants_fun <- function(Source){
  
  df <- ants %>%
    unite(Species2, 
          c("Genus", "Species"),
          sep = " ",
          remove = F) %>%
    mutate(Abundance = as.numeric(Abundance)) %>%
    filter(Source.ID == {{Source}}) %>%
    mutate(Species2 = as.factor(Species2)) %>%
    group_by(Locality.ID) %>%
    complete(Species2, 
             fill = list(Abundance = 0))
  
  return(df)
}

ant_locations <- unique(ants$Source.ID)

ant_list <- lapply(ant_locations, ants_fun)

ant_all <- bind_rows(ant_list) %>%
  ungroup() %>%
  dplyr::select(Abundance) %>%
  rename(COUNT = Abundance) %>%
  mutate(dataset = 'ants')

all_counts <- bind_rows(fish_all, hop_all,
                        bird_all, rodent_all,
                        herps_all, ant_all)

a <- ggplot(all_counts,aes(x = COUNT, fill = dataset)) + 
  geom_histogram(alpha = 0.5) +
  scale_y_sqrt() +
  facet_wrap(~dataset, scales = "free") +
  labs(x = "Number of individuals",
       y = "Number of surveys") +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white"))


a +b
