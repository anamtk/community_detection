# preliminary community metrics stuff
# An Bui
# April 2023

# Load packages ---------------------------------------------------------------

package.list <- c("here", "tidyverse")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

## And loading them

for(i in package.list){library(i, character.only = T)}


# Load matrices -----------------------------------------------------------

# raw presence absence
raw_wide <- readRDS(here("data_outputs",
                         'sbc_fish',
                  "community_matrices",
                  "fish_AQUE1_raw_matrix.RDS")) %>% 
  t() %>% 
  as_tibble(rownames = NA)

raw_long <- raw_wide %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column("year") %>% 
  pivot_longer(cols = 2:52, names_to = "species", values_to = "occupancy") %>% 
  mutate(year_num = case_when(
    year == 2002 ~ 1,
    year == 2003 ~ 2,
    year == 2004 ~ 3, 
    year == 2005 ~ 4,
    year == 2006 ~ 5, 
    year == 2007 ~ 6,
    year == 2008 ~ 7,
    year == 2009 ~ 8,
    year == 2010 ~ 9,
    year == 2011 ~ 10,
    year == 2012 ~ 11,
    year == 2013 ~ 12,
    year == 2014 ~ 13,
    year == 2015 ~ 14,
    year == 2016 ~ 15,
    year == 2017 ~ 16,
    year == 2018 ~ 17,
    year == 2019 ~ 18,
    year == 2020 ~ 19,
    year == 2021 ~ 20,
    year == 2022 ~ 21
  )) %>% 
  select(year, year_num, species, occupancy) %>% 
  arrange(year, species) %>% 
  mutate(year = as.numeric(year))

metadata <- raw_long %>% 
  select(year) %>% 
  unique()

# z-matrix
z <- readRDS(here("data_outputs",
                  'sbc_fish',
                    "community_matrices",
                    "fish_AQUE1_z_matrices.RDS")) %>% 
  mutate(year = as.numeric(year))

wide_listcol <- raw_wide %>% 
  rownames_to_column("year") %>% 
  # some random identifiable number
  mutate(iteration = 2897348) %>% 
  select(iteration, year, ANDA:TSEM) 

z_wide <- z %>% 
  select(iteration, species, year, occupancy) %>% 
  pivot_wider(names_from = "species", values_from = "occupancy") %>% 
  # just to order the iterations in numeric order - i know this doesn't matter
  mutate(iteration = as.numeric(iteration)) %>% 
  arrange(iteration) %>% 
  rbind(wide_listcol)

wide_list <- z_wide %>% 
  nest(data = 2:53) %>% 
  # mutate(data = map(data, ~ select(.x, ANDA)))
  mutate(data = map(data, ~ column_to_rownames(.x, "year"))) %>% 
  mutate(type = case_when(
    iteration == 2897348 ~ "raw",
    TRUE ~ "z"
  ))
