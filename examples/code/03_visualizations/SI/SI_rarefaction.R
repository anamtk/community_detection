

# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 

package.list <- c("here", "tidyverse",
                  "iNEXT", 'patchwork')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

## And loading them

for(i in package.list){library(i, character.only = T)}


# Load data ---------------------------------------------------------------

# Birds -------------------------------------------------------------------

bird <- read.csv(here('examples',
                      'data_output',
                      'bird_fundiv',
                      'tidy_data',
                       'bird_msam_tidy_data.csv'))

bird_yr <- bird %>%
  group_by(RECYEAR, AOUCODE) %>%
  summarise(total = sum(NOBS, na.rm = T)) %>%
  ungroup() %>%
  mutate(total = case_when(total > 0 ~ 1,
                           TRUE ~ 0)) %>%
  ungroup() %>%
  pivot_wider(names_from = RECYEAR,
              values_from = total) %>%
  column_to_rownames(var = 'AOUCODE') %>%
  rowwise() %>%
  mutate(rowsum = sum(c_across(where(is.numeric))))

bird_yr2 <- bird_yr %>%
  dplyr::select(rowsum)

bird_yr2 <- rbind(29, bird_yr2) 

bird_yr2 <- as.vector(bird_yr2)

t <- seq(1, 50, by=1) 

out.bird <- iNEXT(bird_yr2, q=0, datatype="incidence_freq", size=t) 

ggiNEXT(out.bird, type=1, color.var="Assemblage") +  
  theme_bw(base_size = 18) +  
  theme(legend.position="None")+
  labs(x = "Number of sampling years")

bird_cover <- ggiNEXT(out.bird, type=2, color.var="Assemblage") + 
  ylim(c(0.9,1)) + 
  theme_bw(base_size = 18) +  
  theme(legend.position="None") +
  labs(x = "Number of sampling years")


# Grasshoppers ------------------------------------------------------------

hop <- read.csv(here('examples',
                     'data_output',
                     'grasshopper_stability',
                     'tidy_data',
                       'grasshopper_msam_tidy_data.csv'))

hop_yr <- hop %>%
  group_by(YEAR, SPECIES) %>%
  summarise(total = sum(CNT, na.rm = T)) %>%
  ungroup() %>%
  mutate(total = case_when(total > 0 ~ 1,
                           TRUE ~ 0)) %>%
  ungroup() %>%
  pivot_wider(names_from = YEAR,
              values_from = total) %>%
  column_to_rownames(var = 'SPECIES') %>%
  rowwise() %>%
  mutate(rowsum = sum(c_across(where(is.numeric))))

hop_yr2 <- hop_yr %>%
  dplyr::select(rowsum)

hop_yr2 <- rbind(28, hop_yr2) 

hop_yr2 <- as.vector(hop_yr2)

t <- seq(1, 50, by=1) 

out.hop <- iNEXT(hop_yr2, q=0, datatype="incidence_freq", size=t) 

ggiNEXT(out.hop, type=1, color.var="Assemblage") +  
  theme_bw(base_size = 18) +  
  theme(legend.position="None") +
  labs(x = "Number of sampling years")

hop_cover <- ggiNEXT(out.hop, type=2, color.var="Assemblage") + 
  ylim(c(0.9,1)) + 
  theme_bw(base_size = 18) +  
  theme(legend.position="None") +
  labs(x = "Number of sampling years")


# plants ------------------------------------------------------------------

plant <- read.csv(here('examples',
                       'data_output',
                       'plant_turnover',
                       'tidy_data',
                     'plant_msom_tidy_data.csv'))

plant_yr <- plant %>%
  group_by(EventYear, SpecID) %>%
  summarise(total = sum(presence, na.rm = T)) %>%
  ungroup() %>%
  mutate(total = case_when(total > 0 ~ 1,
                           TRUE ~ 0)) %>%
  ungroup() %>%
  pivot_wider(names_from = EventYear,
              values_from = total) %>%
  column_to_rownames(var = 'SpecID') %>%
  rowwise() %>%
  mutate(rowsum = sum(c_across(where(is.numeric))))

plant_yr2 <- plant_yr %>%
  dplyr::select(rowsum)

plant_yr2 <- rbind(10, plant_yr2) 

plant_yr2 <- as.vector(plant_yr2)

t <- seq(1, 25, by=1) 

out.plant <- iNEXT(plant_yr2, q=0, datatype="incidence_freq", size=t) 

ggiNEXT(out.plant, type=1, color.var="Assemblage") +  
  theme_bw(base_size = 18) +  
  theme(legend.position="None")

plant_cover <- ggiNEXT(out.plant, type=2, color.var="Assemblage") + 
  ylim(c(0.8,1)) + 
  theme_bw(base_size = 18) +  
  theme(legend.position="None") +
  labs(x = "Number of sampling years")


# plot together -----------------------------------------------------------

bird_cover /plant_cover /hop_cover +
  plot_annotation(tag_levels = "A")

ggsave(plot = last_plot(),
       filename = here("examples",
                       "pictures",
                       "original_R",
                       "SI_richness_rarefaction.jpg"),
       height = 8,
       width = 4,
       units = "in")
