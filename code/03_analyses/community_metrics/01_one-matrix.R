# preliminary community metrics stuff
# An Bui
# April 2023

# Load packages ---------------------------------------------------------------

package.list <- c("here", "tidyverse", "vegan", "codyn", "patchwork")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

## And loading them

for(i in package.list){library(i, character.only = T)}


# Load matrices -----------------------------------------------------------

# raw presence absence
raw_wide <- readRDS(here("data_outputs",
                  "community_matrices",
                  "fish_AQUE1_raw_matrix.RDS")) %>% 
  t() %>% 
  as_tibble(rownames = NA)

raw_long <- raw_wide %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column("species") %>% 
  pivot_longer(cols = 2:22, names_to = "year", values_to = "occupancy") %>% 
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
  select(species, year, year_num, occupancy) %>% 
  arrange(year, species) %>% 
  mutate(year = as.numeric(year))

metadata <- raw_long %>% 
  select(year) %>% 
  unique()

# z-matrix
z <- readRDS(here("data_outputs",
                    "community_matrices",
                    "fish_AQUE1_z_matrices.RDS")) %>% 
  mutate(year = as.numeric(year))

# select an iteration
z1_wide <- z %>% 
  filter(iteration == 1) %>% 
  # put into wide format for processing
  select(species, year, occupancy) %>% 
  pivot_wider(names_from = "species", values_from = "occupancy") %>% 
  column_to_rownames("year") 

z5_wide <- z %>% 
  filter(iteration == 5) %>% 
  # put into wide format for processing
  select(species, year, occupancy) %>% 
  pivot_wider(names_from = "species", values_from = "occupancy") %>% 
  column_to_rownames("year") 

z10_wide <- z %>% 
  filter(iteration == 10) %>% 
  # put into wide format for processing
  select(species, year, occupancy) %>% 
  pivot_wider(names_from = "species", values_from = "occupancy") %>% 
  column_to_rownames("year") 

z1_long <- z %>% 
  filter(iteration == 4) %>% 
  select(species, year, year_num, occupancy) %>% 
  mutate(year = as.numeric(year))


# calculate species turnover ----------------------------------------------

# using codyn::turnover()
# Computes species turnover between time periods as the proportion of species either gained or lost relative to the total number of species observed across both time periods.
# total turnover = (spp gained + spp lost)/total spp between two timepoints
# disappearance: only species that disappear

# ⊣ total turnover --------------------------------------------------------

raw_total_turnover <- turnover(df = raw_long, time.var = "year", species.var = "species", abundance.var = "occupancy")

z1_total_turnover <- turnover(df = z1_long, time.var = "year", species.var = "species", abundance.var = "occupancy")

# for the whole z-matrix
z_total_turnover <- turnover(df = z, time.var = "year", species.var = "species", abundance.var = "occupancy", replicate.var = "iteration")

total_turnover_plot <- ggplot(data = z_total_turnover, aes(x = year, y = total)) +
  geom_line(aes(group = iteration), alpha = 0.2, color = "darkgreen") +
  geom_line(data = raw_total_turnover, aes(x = year, y = total), color = "#000000") +
  labs(x = "Year", y = "Total species turnover") +
  theme_bw() +
  theme(legend.position = "none")


# ⊣ disappearance ---------------------------------------------------------

raw_dis_turnover <- turnover(df = raw_long, time.var = "year", species.var = "species", abundance.var = "occupancy", metric = "disappearance")

z_dis_turnover <- turnover(df = z, time.var = "year", species.var = "species", abundance.var = "occupancy", replicate.var = "iteration", metric = "disappearance")

dis_turnover_plot <- ggplot(data = z_dis_turnover, aes(x = year, y = disappearance)) +
  geom_line(aes(group = iteration), alpha = 0.2, color = "cornflowerblue") +
  geom_line(data = raw_dis_turnover, aes(x = year, y = disappearance), color = "#000000") +
  labs(x = "Year", y = "Turnover in disappearing species") +
  theme_bw() +
  theme(legend.position = "none")


# ⊣ appearance ------------------------------------------------------------

raw_app_turnover <- turnover(df = raw_long, time.var = "year", species.var = "species", abundance.var = "occupancy", metric = "appearance")

z_app_turnover <- turnover(df = z, time.var = "year", species.var = "species", abundance.var = "occupancy", replicate.var = "iteration", metric = "appearance")

app_turnover_plot <- ggplot(data = z_app_turnover, aes(x = year, y = appearance)) +
  geom_line(aes(group = iteration), alpha = 0.2, color = "orange") +
  geom_line(data = raw_app_turnover, aes(x = year, y = appearance), color = "#000000") +
  labs(x = "Year", y = "Turnover in gained species") +
  theme_bw() +
  theme(legend.position = "none")


# ⊣ putting everything together -------------------------------------------
 
total_turnover_plot / (dis_turnover_plot + app_turnover_plot)


# calculate rate of change ------------------------------------------------

# calculating raw rate change
raw_ratechange <- rate_change_interval(df = raw_long, time.var = "year", species.var = "species", abundance.var = "occupancy")

# double checking against vegdist
raw_eucdist <- vegdist(raw_wide, method = "euclidean")

# plotting rate change
raw_ratechange_plot <- ggplot(data = raw_ratechange, aes(x = interval, y = distance)) +
  geom_point() +
  scale_y_continuous(limits = c(0.8, 5.2)) +
  geom_smooth(se = FALSE, method = "lm", linewidth = 2) +
  labs(x = "Time interval", y = "Euclidean distance", 
       title = "Raw data") +
  theme_bw() 

z_ratechange <- rate_change_interval(df = z, time.var = "year", species.var = "species", abundance.var = "occupancy", replicate.var = "iteration")

z_ratechange_plot <- ggplot(data = z_ratechange, aes(x = interval, y = distance)) +
  geom_point() +
  scale_y_continuous(limits = c(0.8, 5.2)) +
  geom_smooth(se = FALSE, method = "lm", linewidth = 2) +
  labs(x = "Time interval", y = "Euclidean distance", 
       title = "Z") +
  theme_bw() 

raw_ratechange_plot + z_ratechange_plot

# community composition ---------------------------------------------------


# ⊣ raw matrix ------------------------------------------------------------



raw_jacc <- metaMDS(raw_wide, distance = "jaccard")

raw_sitescores <- scores(raw_jacc, "sites") %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column("year")


# ⊣ z-matrix: iteration 1 -------------------------------------------------

z1_jacc <- metaMDS(z1_wide, distance = "jaccard")

z1_sitescores <- scores(z1_jacc, "sites") %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column("year")


# ⊣ z-matrix: iteration 5 -------------------------------------------------

z5_jacc <- metaMDS(z5_wide, distance = "jaccard")

z5_sitescores <- scores(z5_jacc, "sites") %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column("year")


# ⊣ z-matrix: iteration 10 ------------------------------------------------

z10_jacc <- metaMDS(z10_wide, distance = "jaccard")

z10_sitescores <- scores(z10_jacc, "sites") %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column("year")

ggplot(data = raw_sitescores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(size = 3, shape = 24, fill = "#000000") +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(data = z1_sitescores, aes(x = NMDS1, y = NMDS2), alpha = 0.5, color = "blue") +
  geom_point(data = z5_sitescores, aes(x = NMDS1, y = NMDS2), alpha = 0.5, color = "red") +
  geom_point(data = z10_sitescores, aes(x = NMDS1, y = NMDS2), alpha = 0.5, color = "purple") +
  theme_bw()

# beta diversity ----------------------------------------------------------

# note: adespatial::TBI() can do temporal beta diversity for multiple sites at different time points

# a: shared species
# b: species unique to group B
# c: species unique to group C
# unclear... what abc is here. not sure how they all fit together?
raw_abc <- betadiver(raw_wide)

# whittaker beta diversity
raw_whit <- betadiver(raw_wide, "w")

z1_whit <- betadiver(z1_wide, "w")

# my idea is to pull out the diagonals and plot with pair years (e.g. 2002-2003, 2003-2004) on the x axis with beta diversity on the y-axis



