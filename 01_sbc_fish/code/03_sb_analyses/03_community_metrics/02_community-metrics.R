# preliminary community metrics stuff
# An Bui
# April 2023

# source ------------------------------------------------------------------

package.list <- c("here", "tidyverse", "vegan", "codyn", "patchwork")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

## And loading them

for(i in package.list){library(i, character.only = T)}

source(here("code", 
            'sbc_fish',
            "03_sb_analyses",
            "02_community_metrics",
            "01_generating-list-cols.R"))

# calculate species turnover ----------------------------------------------

# using codyn::turnover()
# Computes species turnover between time periods as the proportion of species either gained or lost relative to the total number of species observed across both time periods.
# total turnover = (spp gained + spp lost)/total spp between two timepoints
# disappearance: only species that disappear

# ⊣ total turnover --------------------------------------------------------

raw_total_turnover <- turnover(df = raw_long, time.var = "year", species.var = "species", abundance.var = "occupancy")

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


# ⊣ plotting together -----------------------------------------------------

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


# ⊣ plotting together -----------------------------------------------------

raw_ratechange_plot + z_ratechange_plot

# community composition ---------------------------------------------------

jacc_list <- wide_list %>% 
  mutate(jacc = map(data, ~ metaMDS(.x, distance = "jaccard"))) %>% 
  mutate(jacc_scores = map(jacc, 
                           ~ scores(.x, "sites") %>% 
                             as_tibble(rownames = NA) %>% 
                             rownames_to_column("year")))

# plot raw first
ggplot(data = jacc_list[[5]][[101]], aes(x = NMDS1, y = NMDS2)) +
  geom_point(size = 3, shape = 24, fill = "#000000") +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  # iteration 10
  geom_point(data = jacc_list[[5]][[10]], aes(x = NMDS1, y = NMDS2), alpha = 0.5, color = "red") +
  # iteration 20
  geom_point(data = jacc_list[[5]][[20]], aes(x = NMDS1, y = NMDS2), alpha = 0.5, color = "orange") +
  # iteration 30
  geom_point(data = jacc_list[[5]][[30]], aes(x = NMDS1, y = NMDS2), alpha = 0.5, color = "yellow") +
  # iteration 40
  geom_point(data = jacc_list[[5]][[40]], aes(x = NMDS1, y = NMDS2), alpha = 0.5, color = "green") +
  # iteration 50
  geom_point(data = jacc_list[[5]][[50]], aes(x = NMDS1, y = NMDS2), alpha = 0.5, color = "blue") +
  # iteration 60
  geom_point(data = jacc_list[[5]][[60]], aes(x = NMDS1, y = NMDS2), alpha = 0.5, color = "purple") +
  # iteration 70
  geom_point(data = jacc_list[[5]][[70]], aes(x = NMDS1, y = NMDS2), alpha = 0.5, color = "brown") +
  # iteration 80
  geom_point(data = jacc_list[[5]][[80]], aes(x = NMDS1, y = NMDS2), alpha = 0.5, color = "cornflowerblue") +
  # iteration 90
  geom_point(data = jacc_list[[5]][[90]], aes(x = NMDS1, y = NMDS2), alpha = 0.5, color = "coral") +
  # iteration 100
  geom_point(data = jacc_list[[5]][[100]], aes(x = NMDS1, y = NMDS2), alpha = 0.5, color = "deeppink3") +
  theme_bw() +
  theme(panel.grid = element_blank())

# beta diversity ----------------------------------------------------------

# note: adespatial::TBI() can do temporal beta diversity for multiple sites at different time points

# a: shared species
# b: species unique to group B
# c: species unique to group C
# unclear... what abc is here. not sure how they all fit together?
raw_abc <- betadiver(raw_wide)

# whittaker beta diversity
betadiv_list <- wide_list %>% 
  mutate(whit = map(data, ~ betadiver(.x, "w")))

# my idea is to pull out the diagonals and plot with pair years (e.g. 2002-2003, 2003-2004) on the x axis with beta diversity on the y-axis


# basic diversity metrics -------------------------------------------------

univ_list <- wide_list %>% 
  # calculate shannon diversity for each iteration and raw
  mutate(shandiv = map(data, 
                       ~ diversity(.x, index = "shannon") %>% 
                         as_tibble(rownames = NA) %>% 
                         rownames_to_column("year") %>% 
                         rename("shandiv" = value))) %>% 
  # calculate mean shannon diversity across all years
  mutate(mean_shandiv = map(shandiv, ~ mean(.x$shandiv))) %>% 
  mutate(se_shandiv = map(shandiv, ~ sd(.x$shandiv)/sqrt(length(.x$shandiv)))) %>% 
  # calculate species richness
  mutate(richness = map(data,
                        ~ specnumber(.x) %>% 
                          as_tibble(rownames = NA) %>% 
                          rownames_to_column("year") %>% 
                          rename("richness" = value))) %>% 
  # calculate mean species richness across all years
  mutate(mean_richness = map(richness, ~ mean(.x$richness))) %>% 
  mutate(se_richness = map(richness, ~ sd(.x$richness)/sqrt(length(.x$richness)))) %>% 
  # calculate species frequency
  mutate(frequency = map(data,
                         ~ specnumber(.x, MARGIN = 2) %>% 
                           as_tibble(rownames = NA) %>% 
                           rownames_to_column("species") %>% 
                           rename("frequency" = value))) 


# ⊣ shannon diversity and richness ----------------------------------------

# pulling means out and putting into a data frame to plot
shandiv_mean <- univ_list %>% 
  select(iteration, type, mean_shandiv, se_shandiv) %>% 
  mutate(mean = unlist(mean_shandiv),
         se = unlist(se_shandiv)) %>% 
  mutate(x = "Shannon diversity")

richness_mean <- univ_list %>% 
  select(iteration, type, mean_richness, se_richness) %>% 
  mutate(mean = unlist(mean_richness),
         se = unlist(se_richness)) %>% 
  mutate(x = "Species richness")

# plotting
shandiv_plot <- ggplot(data = shandiv_mean %>% filter(type == "z"), aes(x = x, y = mean)) +
  geom_violin(fill = "aquamarine4", color = "aquamarine4", alpha = 0.2) +
  stat_summary(geom = "point", fun = "median", color = "aquamarine4", size = 3) +
  geom_point(data = shandiv_mean %>% filter(type == "raw"), aes(x = x, y = mean), size = 3) +
  geom_linerange(data = shandiv_mean %>% filter(type == "raw"), aes(x = x, y = mean, ymin = mean - se, ymax = mean + se)) +
  labs(y = "Mean Shannon diversity") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank())

richness_plot <- ggplot(data = richness_mean %>% filter(type == "z"), aes(x = x, y = mean)) +
  geom_violin(fill = "darkorange4", color = "darkorange4", alpha = 0.2) +
  stat_summary(geom = "point", fun = "median", color = "darkorange4", size = 3) +
  geom_point(data = richness_mean %>% filter(type == "raw"), aes(x = x, y = mean), size = 3) +
  geom_linerange(data = richness_mean %>% filter(type == "raw"), aes(x = x, y = mean, ymin = mean - se, ymax = mean + se)) +
  labs(y = "Mean species richness") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank())


# ⊣ ⊣ putting everything together -----------------------------------------

shandiv_plot + richness_plot


# ⊣ species frequency -----------------------------------------------------

ggplot(data = univ_list[[10]][[101]], aes(x = frequency, y = species)) +
  geom_point(size = 3, shape = 24, fill = "#000000") +
  # iteration 10
  geom_point(data = univ_list[[10]][[10]], aes(x = frequency, y = species), alpha = 0.5, color = "red") +
  # iteration 20
  geom_point(data = univ_list[[10]][[20]], aes(x = frequency, y = species), alpha = 0.5, color = "orange") +
  # iteration 30
  geom_point(data = univ_list[[10]][[30]], aes(x = frequency, y = species), alpha = 0.5, color = "yellow") +
  # iteration 40
  geom_point(data = univ_list[[10]][[40]], aes(x = frequency, y = species), alpha = 0.5, color = "green") +
  # iteration 50
  geom_point(data = univ_list[[10]][[50]], aes(x = frequency, y = species), alpha = 0.5, color = "blue") +
  # iteration 60
  geom_point(data = univ_list[[10]][[60]], aes(x = frequency, y = species), alpha = 0.5, color = "purple") +
  # iteration 70
  geom_point(data = univ_list[[10]][[70]], aes(x = frequency, y = species), alpha = 0.5, color = "brown") +
  # iteration 80
  geom_point(data = univ_list[[10]][[80]], aes(x = frequency, y = species), alpha = 0.5, color = "cornflowerblue") +
  # iteration 90
  geom_point(data = univ_list[[10]][[90]], aes(x = frequency, y = species), alpha = 0.5, color = "coral") +
  # iteration 100
  geom_point(data = univ_list[[10]][[100]], aes(x = frequency, y = species), alpha = 0.5, color = "deeppink3") +
  labs(x = "Frequency", y = "Species") +
  theme_bw() +
  theme(panel.grid = element_blank())









