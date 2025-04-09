#relationships between raw and modeled abundances


# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse",
                  'readxl', 'patchwork',
                  'ggridges')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

theme_set(theme_classic())
theme_update(plot.title.position = "plot",
             panel.grid = element_blank(),
             axis.title = element_text(size = 11),
             axis.text = element_text(size = 10),
             plot.title = element_text(size = 12),
             legend.text = element_text(size = 9),
             legend.key.size = unit(0.25, "cm"),
             legend.key = element_blank(), 
             legend.title=element_blank(),
             legend.background = element_blank())

modeled_col <- "#E88C23"
observed_col <- "#438AA8"
# Data --------------------------------------------------------------------

birds <- read.csv(here('examples',
                       'data_output',
                       'bird_fundiv',
                       'tidy_data',
                       'bird_msam_tidy_data.csv')) %>%
  mutate(AOUCODE = case_when(COMMONNAME == "American Goldfinch" ~ 'AGOL',
                             TRUE ~ AOUCODE))

samples <- readRDS(here('examples',
                        'data_output',
                        'bird_fundiv',
                        'computing_cluster_outputs',
                        'bird_N_samples.RDS'))


grasshoppers <- read.csv(here('examples',
                              'data_output',
                              'grasshopper_stability',
                              'tidy_data',
                              'grasshopper_msam_tidy_data.csv')) %>%
  dplyr::select(yrID, siteID, speciesID, CNT, rep)
# 
# #N [speceis, site, year]
grasshopper_modeled <- readRDS(here('examples',
                                    'data_output',
                                    'grasshopper_stability',
                                    'computing_cluster_outputs',
                                    'sev_N_samples.RDS'))

plants <- read.csv(here('examples',
                        'data_output',
                        'plant_turnover',
                        'tidy_data',
                        'plant_msom_tidy_data.csv'))

plant_model <- readRDS(here('examples',
                            'data_output',
                            'plant_turnover',
                            'computing_cluster_outputs',
                            'plant_msom_z_sum.RDS'))
# Prep datasets -----------------------------------------------------------

# Prep modeled N dataset --------------------------------------------------

#n for some species is abnormally large,
#but how to get it down to a reasonable size - 
#maybe sample from the values of that species 
#across other years???

ndf <- as.data.frame.table(samples$N) %>%
  mutate(Iteration = as.numeric(as.factor(Var1)),
         SpecID = as.numeric(as.factor(Var2)),
         TransID = as.numeric(as.factor(Var3)),
         yrID = as.numeric(as.factor(Var4))) %>%
  rename(N = Freq) %>%
  dplyr::select(Iteration, SpecID, TransID, yrID, N)


# Set extremely high values to be more realistic --------------------------

#get the max raw by species plus sd
#when non-observed removed
#set values > than that in model to:
## rnorm(max+sd, sd)
bird_stats <- birds %>%
  filter(NOBS > 0) %>%
  group_by(SpecID) %>%
  summarise(max = max(NOBS, na.rm = T),
            var = var(NOBS, na.rm = T),
            sd = sd(NOBS, na.rm = T)) %>%
  mutate(sd = case_when(is.na(sd) ~ max,
                        TRUE ~ sd),
         var = case_when(is.na(var) ~ max,
                         TRUE ~ var))

ndf2 <- ndf %>%
  left_join(bird_stats, by = "SpecID") %>%
  mutate(N = case_when(N > (max+var) ~ rnorm(1, (max+sd), sd),
                       TRUE ~ N)
  )

ndf3 <- ndf2 %>%
  group_by(SpecID, TransID, yrID) %>%
  summarise(mean_N = mean(N, na.rm = T)) %>%
  ungroup()

# Combine and plot --------------------------------------------------------

bird_df <- birds %>%
  dplyr::select(TransID, yrID, SpecID, NOBS) %>%
  left_join(ndf3, by = c("TransID", "yrID", "SpecID")) %>%
  pivot_longer(NOBS:mean_N,
               names_to = "type",
               values_to = "N")
# 
# ggplot(bird_df, aes(x = N, fill = type)) +
#   geom_density() +
#   scale_x_sqrt() +
#   scale_y_sqrt()

a <- ggplot(bird_df, aes(x = N, fill = type)) +
  geom_histogram(position = "dodge") +
  scale_x_sqrt() +
  scale_y_sqrt() +
  scale_fill_manual(labels = c("mean_N" = "modeled",
                               "NOBS" = "observed"),
                      name = "",
                      values = c(modeled_col, 
                                 observed_col)) +
  labs(x = "Number of individuals\n(per species)",
       y = "Number of surveys")



# Grasshoppers ------------------------------------------------------------

gndf <- as.data.frame.table(grasshopper_modeled$N) %>%
  mutate(Iteration = as.numeric(as.factor(Var1)),
         speciesID = as.numeric(as.factor(Var2)),
         siteID = as.numeric(as.factor(Var3)),
         yrID = as.numeric(as.factor(Var4))) %>%
  rename(N = Freq) %>%
  dplyr::select(Iteration, speciesID, siteID, yrID, N)

hopper_sum <- grasshoppers %>%
  filter(CNT > 0) %>%
  group_by(speciesID) %>%
  summarise(max = max(CNT, na.rm = T),
            sd = sd(CNT, na.rm = T),
            var = var(CNT, na.rm = T)) %>%
  mutate(sd = case_when(is.na(sd) ~ max,
                        TRUE ~ sd),
         var = case_when(is.na(var) ~ max,
                         TRUE ~ var))

gndf2 <- gndf %>%
  left_join(hopper_sum, by = "speciesID")%>%
  mutate(N = case_when(N > (max+var) ~ rnorm(1, (max+sd), sd),
                       TRUE ~ N)) 

gndf3 <- gndf2 %>%
  group_by(siteID, yrID, speciesID) %>%
  summarise(mean_N = mean(N, na.rm = T)) %>%
  ungroup()

hop_df <- grasshoppers %>%
  dplyr::select(siteID, yrID, speciesID, CNT) %>%
  left_join(gndf3, by = c('siteID', 'yrID', 'speciesID')) %>%
  pivot_longer(CNT:mean_N,
               names_to = "type",
               values_to = "N")
# 
# ggplot(hop_df, aes(x = N, fill = type)) +
#   geom_density(alpha = 0.5) +
#   scale_x_sqrt() +
#   scale_y_sqrt() +
#   scale_fill_manual(labels = c("mean_N" = "modeled",
#                              "CNT" = "observed"),
#                   name = "",
#                   values = c(observed_col,
#                              modeled_col
#                   )) +
#   labs(x = "Number of individuals\n(per species)",
#        y = "Number of surveys")

b <- ggplot(hop_df, aes(x = N, fill = type)) +
  geom_histogram(position = "dodge") +
  scale_x_sqrt() +
  scale_y_sqrt() +
  scale_fill_manual(labels = c("mean_N" = "modeled",
                               "CNT" = "observed"),
                    name = "",
                    values = c(observed_col,
                               modeled_col
                               )) +
  labs(x = "Number of individuals\n(per species)",
       y = "Number of surveys") +
  theme(legend.position = "none")




# Plants ------------------------------------------------------------------

#z[k,t,i]


plant_model2 <- as.data.frame(plant_model$quantiles) %>%
  rownames_to_column(var = "par") %>%
  filter(par != "deviance") %>%
  mutate(par = str_sub(par, 3, (nchar(par)-1))) %>%
  separate(par,
           into = c("species", "year", "site"),
           sep = ",") %>%
  dplyr::select(species, year, site, `50%`) %>%
  rename(presence = `50%`) 

plant_model2 %>%
  ungroup() %>%
  distinct(year, site) %>%
  tally()

plant_model3 <- plant_model2 %>%
  filter(presence > 0) %>%
  group_by(species) %>%
  summarise(total = n(),
            prop = total/360) %>%
  mutate(type = "modeled")

plants %>%
  ungroup() %>%
  distinct(yrID, Plot, Transect, quadnum) %>%
  tally()

plant2 <- plants %>%
  filter(presence > 0) %>%
  group_by(SpecID)%>%
  summarise(total = n(),
            prop = total/360) %>%
  mutate(type = "observed")

all_plants <- plant_model3 %>%
  bind_rows(plant2)

c <- ggplot(all_plants, aes(x = prop, fill = type)) +
  geom_histogram(position = "dodge")+
  scale_fill_manual(name = "",
                    values = c(modeled_col, 
                               observed_col)) +
  labs(x = "Frequency of observation \n(proportion of surveys)",
       y = "Number of species") +
  theme(legend.position = "none")

# ggplot(all_plants, aes(x = prop, fill = type)) +
#   geom_density(alpha = 0.5)+
#   scale_fill_manual(name = "",
#                     values = c(modeled_col, 
#                                observed_col)) +
#   labs(x = "Presence/absence status",
#        y = "Number of observations")


# Pull together -----------------------------------------------------------

a +  c +b +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a",
                  tag_suffix = ")")

ggsave(here('examples',
            'pictures',
            'original_R',
            'Fig2_what_does_detection_correction_do.pdf'),
       height = 3, 
       width = 10,
       units = "in")

