#relationships between raw and modeled abundances


# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse",
                  'readxl', 'patchwork')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

theme_set(theme_bw())

# Data --------------------------------------------------------------------

birds <- read.csv(here('02_konza_birds',
                       'data_outputs',
                       'MSAM',
                       'knz_tidy_data_for_model.csv')) %>%
  mutate(AOUCODE = case_when(COMMONNAME == "American Goldfinch" ~ 'AGOL',
                             TRUE ~ AOUCODE))

samples <- readRDS(here('02_konza_birds',
                        'monsoon',
                        'MSAM',
                        'outputs',
                        'bird_N_samples.RDS'))




# Prep datasets -----------------------------------------------------------

#get average species-level abundances
bird_abund <- birds %>%
  group_by(SpecID) %>%
  summarise(mean = mean(NOBS),
            sd = sd(NOBS),
            total = n(),
            se = sd/sqrt(total))%>%
  mutate(type = "observed")

ndf <- as.data.frame.table(samples$N) %>%
  mutate(Iteration = as.numeric(as.factor(Var1)),
         SpecID = as.numeric(as.factor(Var2)),
         TransID = as.numeric(as.factor(Var3)),
         yrID = as.numeric(as.factor(Var4))) %>%
  rename(N = Freq) %>%
  dplyr::select(Iteration, SpecID, TransID, yrID, N)

samples_abund <- ndf %>%
  group_by(SpecID) %>%
  summarise(mean = mean(N, na.rm = T),
            sd = sd(N, na.rm = T),
            total = n(),
            se = sd/sqrt(total)) %>%
  mutate(type = "modeled")


# Combine and plot --------------------------------------------------------

obs <- bind_rows(bird_abund, samples_abund)

ggplot(obs, aes(x = mean, fill = type)) +
  geom_density(alpha = 0.5) +
  scale_x_sqrt()

ggplot(obs, aes(x = mean, fill = type)) +
  geom_histogram(position = "dodge") +
  scale_x_sqrt()

# plot raw ----------------------------------------------------------------

bird_site <- birds %>%
  dplyr::select(NOBS, TransID, yrID, SpecID) %>%
  mutate(type = "observed") %>%
  rename(COUNT = NOBS)

samp_site <- ndf %>%
  group_by(SpecID, TransID, yrID) %>%
  summarise(COUNT = mean(N, na.rm = T)) %>%
  mutate(type = "modeled")

rawdf <- bind_rows(bird_site, samp_site)

ggplot(rawdf, aes(x = COUNT, fill = type)) +
  geom_density() +
  scale_x_log10()

ggplot(rawdf, aes(x = COUNT, fill = type)) +
  geom_density() +
  scale_x_sqrt()

ggplot(rawdf, aes(x = COUNT, fill = type)) +
  geom_histogram(position = "dodge") +
  scale_x_sqrt() +
  scale_y_sqrt()+
  labs(x = "Number of individuals",
       y = "Number of observations")

